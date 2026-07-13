.. _operator-algebra:

=============================
Operator Algebra Architecture
=============================

ORPHEUS' transport, eigenvalue, and Krylov solvers all act on a flux
distribution by composing a small set of linear operators — the
invertible loss operator :math:`A = L + C` (streaming :math:`L =
\Omega\cdot\nabla` plus collision :math:`C = \Sigma_t`), the scattering
operator :math:`S`, the fission operator :math:`F`. The
:mod:`orpheus.numerics.operator` module installs these as a uniform
*matrix-free* algebra so that the eigenvalue, fixed-source, and
preconditioned Krylov code can consume any solver method (SN / MoC /
CP / diffusion) without knowing which transport discretisation lives
underneath. This page is the **Phase 0** stub for that algebra; the
full unifying narrative lands progressively under the Wave H refactor
(Issue 18).

.. contents::
   :local:
   :depth: 2


Key Facts
=========

- **The realized boundary law** :math:`B` **is a first-class sibling
  operator** (Wave O step O.4a.2, Issue #208, 2026-06-03): the
  canonical SN transport algebra is now
  :math:`(L_{\rm full} + C - S - F - B)\,\psi = q` on the direct-sum
  transport state :math:`V = V_{\rm bulk} \oplus V_{\rm inflow} \oplus
  V_{\rm outflow}`. The boundary reflection is **no longer re-applied
  inside the streaming sweep** (the deleted "keystone"); it is
  delivered as the off-diagonal :math:`-B` source term and the outer
  Krylov / SI loop drives the boundary consistency residual
  :math:`\psi.\text{inflow} - B\,\psi.\text{outflow} - q.\text{inflow}
  \to 0`. See :ref:`bc-extraction`.

- **Every operator's** ``.apply`` **output boundary is a**
  :class:`~orpheus.transport.source_sinks.angular_boundary_source_sink.AngularBoundarySourceSink`
  (Wave O step B.5.2, Issue #208, ``6ef5063``, 2026-06-03), completing
  the boundary half of the operator-output "dimensional-sin" carve (the
  bulk half — ``.apply.bulk`` →
  :class:`~orpheus.transport.source_sinks.angular_source_sink.AngularSourceSink`
  — landed in ``f400743``). The governing principle: an operator output
  is :math:`A\psi` — a **source/sink**, NOT a residual; the residual
  arises ONLY from an explicit :meth:`~orpheus.transport.residuals.angular_boundary_residual.AngularBoundaryResidual.from_balance`
  of the output against a source. The completed boundary role grid
  mirrors the bulk:
  :math:`\texttt{.apply} \to \texttt{AngularBoundarySourceSink}`,
  :math:`\texttt{.solve} \to \texttt{AngularBoundaryFlux}`,
  :math:`\texttt{from\_balance} \to \texttt{AngularBoundaryResidual}`. See
  :ref:`bc-extraction-operator-output-typing`.

- **Flux states form an affine space; the iterate increment is a typed
  displacement, and** ``flux + flux`` **is a** :class:`TypeError`
  (Wave O step O.2 close-out, Issue #208/#201, 2026-06-08). The flux
  states :class:`~orpheus.transport.fields.angular_flux.AngularFlux` /
  :class:`~orpheus.transport.fields.scalar_flux.ScalarFlux` /
  :class:`~orpheus.transport.fields.harmonic_moment_flux.HarmonicMomentFlux`
  / :class:`~orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux` are an
  **affine space** :math:`A` over a distinct **difference vector space**
  :math:`V`. The source-iteration increment
  :math:`\Delta\psi = \psi^{(i)} \ominus \psi^{(i-1)}` is an element of
  :math:`V` — a :class:`~orpheus.transport.displacements._displacement.Displacement`,
  NOT a state — so ``flux ⊖ flux`` mints a displacement, ``flux ⊕
  displacement`` is the torsor update, and ``flux + flux`` raises (an
  affine space has no origin; the #201 dimensional gate is now a
  **type** consequence). The displacement is the only object that knows
  "previous"/"step", so it carries the convergence diagnostics a flux
  state cannot (``contraction_ratio`` :math:`\approx c`,
  ``true_error_estimate`` = the :math:`c\to 1` false-convergence fix).
  The equation residual :math:`r = (L+C-S-B)\psi - q` is now typed via
  :func:`~orpheus.sn.solver.evaluate_residual` (the box-7 consumer of
  the previously-unconsumed ``from_balance`` mint). See
  :ref:`affine-typed-field-algebra`.

- **The carriers form a** :math:`(\text{Representation} \times
  \text{Role})` **double category, and the operator algebra traverses
  it** (Frame-projection campaign P4.5, branch
  ``refactor/operator-inverse-algebra``, 2026-06-25). A carrier is a cell
  :math:`(\text{Representation}, \text{Role})`: **Representation**
  :math:`\in \{\text{Angular}, \text{Moment}, \text{Scalar},
  \text{Trace}\}` sets the array shape and carries the change-of-basis
  (the Frame); **Role** :math:`\in \{\text{Flux}, \text{Source},
  \text{Residual}, \text{Displacement}\}` sets the arithmetic interface
  (the #208 affine torsor). The **horizontal** 1-morphisms are the
  representation-changing frame faces :math:`M`/:math:`R` (role-generic —
  a base change that fixes the fiber); the **vertical** 1-morphisms are
  the role-changing cross sections :math:`C`/:math:`\Lambda`/:math:`F`
  (representation-generic — the role change *is* the cross-section
  physics); **scattering** :math:`S = \tfrac{1}{W}(R\circ\Lambda\circ M)
  = \texttt{frame.conjugate}(\Lambda)` is the **2-cell**, the vertical
  :math:`\Lambda` conjugated by the horizontal adjoint pair, and the
  bit-identical windowed-vs-full crosscheck is its **interchange-law
  coherence witness**. **A grid cell IS an operator's** ``(Domain,
  Codomain)``:
  :class:`~orpheus.numerics.operator.LinearOperator` ``[Domain,
  Codomain]`` is the typed grid traversal — the parametrization belongs
  on the *operator*, not the carrier, because a fully-typed
  ``Carrier[Representation, Role]`` is **structurally impossible** (Role
  changes ``__add__`` ⟹ Role must be a class; Representation changes
  shape ⟹ Representation must be a class; a parameterized carrier would
  break the runtime units gate via generic erasure). The flat
  multiple-inheritance leaves ``AngularFlux(FluxRole, AngularField)`` are
  therefore the **unique principled normal form**, not a compromise. See
  :ref:`carrier-grid-double-category` and
  :ref:`carrier-grid-flat-leaf-normal-form`.

- **The interior cell-face angular fluxes are a 1-cochain**
  :math:`C^1_{\rm int}` (Wave O step #205 Phase 5, Issue #208,
  2026-06-04): the 2-D wavefront sweep and matvec no
  longer carry raw ephemeral ``psi_x`` / ``psi_y`` numpy arrays. The
  interior 1-cochain :math:`C^1_{\rm int}` and the boundary trace
  :class:`~orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux`
  (the boundary 1-cochain :math:`C^1_\partial`) **biproduct-decompose
  the full face cochain** :math:`C^1 = C^1_{\rm int} \oplus
  C^1_\partial` — the :math:`V_{\rm bulk} \oplus V_{\rm boundary}`
  shape of :eq:`bc-extraction-direct-sum-state` one locus down, at the
  *face* level. The seed/absorb the sweep applies are the typed trace
  operators :math:`\iota_*` / :math:`\iota^*`, with the "absorption =
  identity" fact the provable biproduct law :math:`\iota^* \circ
  \iota_* = \mathrm{id}`. The dedicated ``WavefrontFlux`` carrier was
  **retired at S6.4(f)** (#222): the cochain now lives in the rolling
  front (``_MovingFrontier``) and the full-cochain oracle history
  (``_octant_face_cochain``). See :ref:`wavefront-flux-cochain` for the
  succession.

- **The 2-D Cartesian SI iterate lives in moment space** (Wave O step
  #205 **Phase 5a**, ``93807aa`` / ``b97d4f9`` / ``13ca001``,
  2026-06-07): the within-group source-iteration fixed point
  :math:`\psi_{k+1} = (L{+}C)^{-1}(S\psi_k + B\psi_k + q)` consumes
  :math:`\psi` only through its flux moments :math:`\phi_\ell^m =
  (M\psi)_\ell^m`, so the *persistent* iterate is held as the moment
  tensor :class:`~orpheus.transport.fields.harmonic_moment_flux.HarmonicMomentFlux`
  (:math:`N \to (L{+}1)(2L{+}1)`, measured **18.3×** shrink at
  :math:`N=110`, :math:`L=1`) rather than the full per-ordinate
  :class:`~orpheus.transport.fields.angular_flux.AngularFlux`. The
  source is **bit-identical** (the moment arm of
  :meth:`ScatteringOperator.apply <orpheus.transport.operators.scattering.ScatteringOperator.apply>`
  shares the :math:`R\,\Lambda` reconstruction with the full-angular
  arm); only the SI convergence test moves to the moment :math:`L^2`
  (principled-equivalence). 2-D Cartesian only (curvilinear's
  Morel–Montry Carlson seed reads the per-ordinate iterate; Krylov
  iterates the full bulk). Interior-bulk only — the trace
  :math:`C^1_\partial` stays un-reduced. See :ref:`sn-angular-windowing`.

- **2-D Cartesian eigenvalue problems solve via BOTH inner solvers**
  (Wave O "2-D SI Phase A", Issue #208, 2026-06-04): the
  source-iteration inner
  (:meth:`SNSolver._solve_source_iteration <orpheus.sn.solver.SNSolver._solve_source_iteration>`,
  the :func:`~orpheus.sn.solver.solve_sn` default for *every* geometry)
  AND the Krylov inner
  (:meth:`SNSolver._solve_krylov <orpheus.sn.solver.SNSolver._solve_krylov>`).
  The SI inner is the **geometry-agnostic structural twin** of Krylov:
  identical composite RHS, identical loss decomposition (the invertible
  resolvent :math:`L + C` plus the two lagged coupling gains — the bulk
  scattering :math:`S` and the trace boundary reflection :math:`B` —
  delivered to the **variadic** driver per :ref:`bc-extraction-variadic-driver`;
  zero within-group fission), identical angular reduction — differing
  **only** in the iteration driver. The reflective coupling rides the **bare** 2-D
  sweep via the sibling :math:`-B` on the natively four-face
  (``xmin`` / ``xmax`` / ``ymin`` / ``ymax``)
  :class:`~orpheus.sn.operators.boundary.SNBoundaryOperator` — the same
  operator the 2-D Krylov path uses. The legacy "B1'' face block"
  (never a code symbol; a 1-D boundary-closure *name*) is retired. See
  :ref:`bc-extraction-2d-si-krylov-twin`.

- The transport eigenvalue and fixed-source problems share a single
  operator algebra (Trefethen & Bau 1997, §3.2):

  .. (vv-status rationale) Phase 0 stub label for the transport
     fixed-source operator-algebra form. The verifiable claim — that
     a Phase-1 SN/MoC/CP solver assembles its operator triple
     (A, S, F) and lets ``op = A - S - F`` agree with the legacy
     fixed-source path bit-for-bit — is queued for Issues 9-13.
  .. vv-status: operator-fixed-source documented

  .. math::
     :label: operator-fixed-source

     (A - S - F)\,\psi \;=\; q
     \qquad \text{(fixed source)}

  .. (vv-status rationale) Phase 0 stub label for the transport
     eigenvalue operator-algebra form. Verified by Issue 14's unified
     iteration driver consuming (A - S, F) triples, and by the
     existing eigenvalue test suites once they are migrated to the
     new contract.
  .. vv-status: operator-eigenvalue documented

  .. math::
     :label: operator-eigenvalue

     (A - S)\,\psi \;=\; \tfrac{1}{k}\,F\,\psi
     \qquad \text{(eigenvalue)}

  Both are built from operator addition, scalar multiplication, and
  composition (``+``, ``-``, ``*``, ``@``) acting on
  :class:`~orpheus.numerics.operator.LinearOperator` instances.

- **The eigenvalue problem is layered into four tiers — leaves,
  posing, resolvent, algorithm** (2026-06-05, ``650032e`` /
  ``7603c8e``). The canonical **standard form** is the generalized
  eigenproblem :math:`A_{\rm loss}\,\psi = \lambda\,M\,\psi`, whose
  power-method realization is the dominant eigenpair of the
  **resolvent** :math:`A_{\rm loss}^{-1} M`. The **k-eigenvalue** row
  is :math:`A_{\rm loss} = L+C-S-B`, :math:`M = F`, :math:`k = \mu`;
  the :math:`\alpha`-eigenvalue, adjoint, and transient rows are
  **documented future seams**.
  :func:`~orpheus.numerics.eigenvalue.power_iteration` is the
  **canonical Layer-4 algorithm** (NOT deprecated — it is the *more
  general* layer, binding the resolvent late through the opaque
  :class:`~orpheus.numerics.eigenvalue.EigenvalueSolver` Protocol so it
  admits monolithic-matrix resolvents that have no :math:`(A,S,F)`
  triple);
  :class:`~orpheus.numerics.iteration.KEigenvalue` is the
  operator-triple realization that *delegates* its loop to it (one
  loop in the codebase, Cardinal Rule 2). See
  :ref:`eigenvalue-posing`.

- **The Hilbert adjoint** ``op.H`` **is the metric-correct G-adjoint**
  :math:`A^{\dagger} = G^{-1} A^{\mathsf T} G` (step O.2b R5,
  ``5c06196``), NOT the Euclidean transpose. For the SN composite the
  metric :math:`G` is **block-diagonal** on :math:`V_{\rm
  bulk}\oplus V_{\rm trace}`: bulk :math:`V_{\rm cell}\,w_n` (phase-space
  measure) :math:`\oplus` trace :math:`|\Omega\cdot\hat n_f|\,w_n`
  (partial-current surface measure, pseudo-inverted on the singular
  tangential ordinates). The carrier is
  :class:`~orpheus.numerics.spaces.full_field_space.FullFieldSpace`;
  ``L``/``C``/``S``/``F``/``B`` all carry it (P4.5 W-D) so the
  within-group :class:`~orpheus.numerics.operator.OperatorSum` guard
  VALIDATES the composition, and — because every loss leaf carries the
  composite metric — the adjoint is applied **once at the op level**
  (the ``_AdjointOperator`` wrapper reads :math:`G` off the composite
  domain) and is never metric-blind. Any non-adjointable operand still
  makes its composite non-adjointable, so the recursive
  :attr:`~orpheus.numerics.operator.LinearOperator.is_adjointable`
  reports ``False`` and ``.H`` *raises*
  :class:`~orpheus.numerics.operator.MissingAdjoint` **eagerly at
  construction**, never silently goes Euclidean. See
  :ref:`bc-extraction-g-adjoint`.

- The :class:`~orpheus.numerics.operator.LinearOperator` Protocol
  carries one mandatory method (``apply``) and, per optional axis
  (inverse, adjoint, and — since stencil-assembly 2b — **assembly**,
  :ref:`operator-algebra-assembly-axis`), a **three-layer structural
  surface** (#226 carve P4, Design C): a runtime **predicate**
  (:attr:`~orpheus.numerics.operator.LinearOperator.is_invertible` /
  :attr:`~orpheus.numerics.operator.LinearOperator.is_adjointable` /
  :attr:`~orpheus.numerics.operator.LinearOperator.is_assemblable`), an
  **operator-returning method** (``inverse()`` per-class /
  :attr:`~orpheus.numerics.operator.LinearOperator.H` on the base /
  ``assemble()`` →
  :class:`~orpheus.numerics.assembled_operator.SparseAssembledOperator`),
  and a **realization verb** (``solve`` / ``apply_transpose``) present
  only where a native realization exists — each axis refusing eagerly via
  its own ``TypeError`` sibling
  (:class:`~orpheus.numerics.operator.NotInvertible` /
  :class:`~orpheus.numerics.operator.MissingAdjoint` /
  :class:`~orpheus.numerics.operator.MissingAssembly`). The stringly-typed
  ``capabilities: frozenset[str]`` advertisement it replaced (``CAP_*``
  tags + ``MissingCapability``) is **retired** — the surface itself is
  now the single source of truth (:ref:`capability-set-semantics`).

- **SN array-storage convention** for every operator leaf
  (:class:`~orpheus.sn.operators.streaming.StreamingOperator`, the collision
  multiplier :math:`C = M[\sigma_t]`
  (:class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`),
  :class:`~orpheus.transport.operators.scattering.ScatteringOperator`,
  :class:`~orpheus.transport.operators.fission.FissionOperator`): the
  ``apply(psi) -> psi'`` contract consumes and returns
  ``psi.shape == (N, ng, nx, ny)`` for angular flux,
  ``phi.shape == (ng, nx, ny)`` for scalar flux.  The canonical
  statement with derivation and migration history lives at
  :ref:`theory-sn-index-convention`.

- The **operator surface itself** is the single source of truth for
  what an operator can do — a method exists (or is refused) exactly
  where the ability does, with no parallel string registry that could
  drift. Composition primitives compute their own
  :attr:`~orpheus.numerics.operator.LinearOperator.is_invertible` /
  :attr:`~orpheus.numerics.operator.LinearOperator.is_adjointable`
  recursively from the operands; mismatches fail at composition time,
  NEVER mid-iteration — an ``apply`` mismatch raises ``TypeError``
  eagerly, a non-adjointable ``.H`` raises
  :class:`~orpheus.numerics.operator.MissingAdjoint` at construction,
  and a value-dependent ``inverse()`` raises
  :class:`~orpheus.numerics.operator.NotInvertible` before any inverse
  object exists.

- **The curvilinear ψ½ ray is System B of a 2×2 coupled block operator**
  (coupled-block campaign, GH #280/#282, branch
  ``refactor/sn-walk-unification``, 2026-07-12). The within-group
  augmented S\ :sub:`N` problem is posed as
  :math:`\bigl[\begin{smallmatrix} A_{AA} & A_{AB} \\ A_{BA} & A_{BB}
  \end{smallmatrix}\bigr]\bigl[\begin{smallmatrix} \psi_A \\ \psi_B
  \end{smallmatrix}\bigr] = \bigl[\begin{smallmatrix} q_A \\ q_B
  \end{smallmatrix}\bigr]` over **System A** (the transport
  :class:`~orpheus.transport.full_field.FullField`) and **System B** (the
  ψ½ radial-characteristic ray, its own
  :class:`~orpheus.transport.radial_characteristic.RadialCharacteristicField`
  closed by a two-point BVP). The four blocks are named operators in
  :mod:`orpheus.sn.operators.radial_characteristic` — the seed
  :math:`A_{AB}` (σ-independent), the emission
  :math:`A_{BA} = \mathrm{Fold}\circ K_{\rm iso}\circ\int d\mu`, the
  radial march :math:`A_{BB}` (its direct Carlson solve IS the inverse) —
  assembled by the N-general
  :class:`~orpheus.numerics.coupled_system.CoupledOperator` machinery. The
  one production spelling is
  :func:`~orpheus.sn.coupled_system.build_within_group_system`, which
  returns the :class:`~orpheus.sn.coupled_system.WithinGroupSystem` record
  carrying the named splitting :math:`A = M - N`: the resolvent
  :math:`M = \bigl[\begin{smallmatrix} L+C & \text{Seeding} \\ \mathbf 0
  & A_{BB}\end{smallmatrix}\bigr]` solves block-triangular (System B
  first), the emission gain rides :math:`N` (lagged). Presence is
  **structural** (R12a — System B exists iff the mesh carries a ray; a
  mismatched composite is unconstructable). The stop is the ρ-honest
  free-identity residual with a driver-level lag-death
  :class:`~orpheus.sn.solver.ConvergenceCertificateError`. See
  :ref:`coupled-block-operator`.


Definitions
===========

The three primitive actions on a flux vector :math:`\psi`:

.. (vv-status rationale) Phase 0 stub label for the ``apply`` primitive
   action. Verified at the protocol level by
   ``tests/numerics/test_operator.py`` (foundation-tagged software
   invariants); a per-solver Phase-1 test will check that each
   solver's ``L.apply(x)`` matches its legacy path bit-for-bit.
.. vv-status: operator-apply documented

.. math::
   :label: operator-apply

   \texttt{apply} \;:\; x \;\mapsto\; L\,x

.. (vv-status rationale) Phase 0 stub label for the ``solve``
   primitive action — the algorithmic dual of ``apply``, NOT the
   matrix inverse. Verified at the protocol level by
   ``tests/numerics/test_operator.py``; per-solver verification is
   queued for the BiCGSTAB-consumer migration (Issue 15).
.. vv-status: operator-solve documented

.. math::
   :label: operator-solve

   \texttt{solve} \;:\; b \;\mapsto\; L^{-1}\,b
   \quad \text{(the algorithmic dual of } \texttt{apply}\text{)}

.. (vv-status rationale) Phase 0 stub label for the
   ``apply_transpose`` primitive action. Verified at the protocol
   level by ``tests/numerics/test_operator.py``; per-solver adjoint
   sensitivity tests are queued for the sensitivity track (Issue 17).
.. vv-status: operator-apply-transpose documented

.. math::
   :label: operator-apply-transpose

   \texttt{apply\_transpose} \;:\; x \;\mapsto\; L^{T}\,x

The dual relationship in :eq:`operator-solve` is **algorithmic**, not
matrix-theoretic: ``solve(L, b)`` returns whatever vector the
operator's solve algorithm produces such that
``apply(solve(b)) ≈ b`` to working precision. For a sparse direct
factorisation this is the matrix inverse; for an iterative
preconditioned solver it is an approximate inverse; for a
matrix-vector product wrapping BiCGSTAB it is the Krylov approximate
solution. The protocol does NOT require ``apply ∘ solve = I`` exactly
— only that ``solve`` realizes a meaningful approximation. Operators
without an efficient inverse action simply report
:attr:`~orpheus.numerics.operator.LinearOperator.is_invertible`
``= False`` and carry no ``solve`` verb, and downstream code declines
them at composition time (:ref:`capability-set-semantics`).


.. _heteromorphic-apply-typing:

Typing the heteromorphic ``apply`` — Pattern M
----------------------------------------------

The :eq:`operator-apply` contract is nominally an **endomorphism**
``apply(x: V) -> V`` — flux in, flux of the same type out. That is the
honest signature for the streaming, collision, and boundary leaves. But
two operators are **heteromorphic multi-carrier**: their ``apply``
maps *each input carrier to a distinct output carrier*. The fission and
scattering operators (:ref:`integral-kernel-category`) dispatch on the
*input* carrier type via :func:`functools.singledispatchmethod` and map,
for example,
:class:`~orpheus.transport.fields.scalar_flux.ScalarFlux` to
:class:`~orpheus.transport.source_sinks.ScalarSourceSink`, the timeless
composite :class:`~orpheus.transport.full_field.FullField` to a
(source-role) :class:`~orpheus.transport.full_field.FullField`, and
(scattering)
:class:`~orpheus.transport.fields.harmonic_moment_flux.HarmonicMomentFlux`
to :class:`~orpheus.transport.source_sinks.AngularSourceSink`. The
composite arm registers on (and the ``@overload`` surface names) the
**timeless** :class:`~orpheus.transport.full_field.FullField` — the
history-blind operator carrier (P4.5 W-C/W-F); a driver's history-bearing
:class:`~orpheus.transport.timed_full_field.TimedFullField` iterate reaches
it via MRO (it *is* a ``FullField``). The output
type is **not** ``V`` — it is a function of the input carrier (the §B.5.2
truth that an operator output is a *source/sink*, not a flux). Naming
these maps honestly to the type checker is the first use of
:func:`typing.overload` anywhere in ``orpheus/`` (#257 S8c), and it
required a small idiom worth documenting because it will recur.

Why the obvious spellings fail
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Python has no native multiple dispatch, and the two stdlib tools each
fall short alone:

* **A raising endomorphism base poisons the type.** The natural
  ``@singledispatchmethod`` shape has a *base* method that raises
  :class:`TypeError` for unregistered carriers. If that base is named
  ``apply`` and inherits the mixin's nominal ``apply(x: V) -> V``, the
  raising body makes pyright infer ``singledispatchmethod[NoReturn]`` —
  so every caller of ``apply`` statically sees ``NoReturn`` (a poison
  type that contaminates downstream inference), and every
  ``@apply.register`` arm errors against the inherited endomorphism
  signature.

* **An** :func:`~typing.overload` **stub cannot carry the body.**
  ``@overload`` is a type-checker fiction erased at runtime: an overload
  signature is *only* a signature, so it can never hold the
  source-assembly math (≈150 lines for scattering). And
  ``@singledispatchmethod`` *does* carry per-carrier bodies with
  automatic routing — its only failing is the homomorphic
  ``singledispatchmethod[_T]`` typing.

Pattern M — the overload surface over the renamed dispatcher
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The fix (**Pattern M**) keeps ``@singledispatchmethod`` for its routing
and bodies, and bolts a typed surface on top:

1. **Rename the dispatcher** ``apply`` → ``_apply_impl`` with a
   ``-> Any`` base, so each ``@_apply_impl.register`` arm is a bodied,
   *real-typed* function (e.g. ``def _(self, phi: ScalarFlux) ->
   ScalarSourceSink``) that ``.register`` accepts at its natural
   indentation.

2. **Add the typed surface only for the type checker:**

   .. code-block:: python

      if TYPE_CHECKING:
          @overload
          def apply(self, phi: ScalarFlux, /) -> ScalarSourceSink: ...
          @overload
          def apply(self, psi: FullField, /) -> FullField: ...
          # ... one overload per carrier ...
          def apply(self, x: Any, /) -> Any: ...
      else:
          apply = _apply_impl

   At runtime the ``else`` branch makes the public ``apply`` the **same
   object** as the dispatcher — ``Type.__dict__['apply'] is
   Type.__dict__['_apply_impl']`` is ``True``. The public ``apply`` *is*
   the ``singledispatchmethod``, merely aliased: **runtime is
   byte-identical** to the untyped version, while pyright sees the
   per-carrier overloads and reports the exact output type for each
   input.

Why Pattern M over the alternatives
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Pattern M was chosen over a ``TYPE_CHECKING`` / ``else`` *split* of the
whole method (call it Pattern C — defining a fully-typed ``apply`` under
``if TYPE_CHECKING`` and the dispatcher under ``else``). The decision is
the master-standard rule ordering: "code reveals intention" (Beck rule 2)
outranks "fewest elements" (rule 4). Pattern C buries the bulk of the
source-assembly math inside an ``else:`` block (a visual demotion of the
operator's actual work); Pattern M extends the mixin's *existing*
``if TYPE_CHECKING: def apply(self, x: V) -> V`` idiom (an
internally-consistent precedent), keeping the dispatch arms and their
bodies at natural indentation where a reader expects to find them.

The deeper *spelling* question — Pattern M (``@singledispatchmethod`` +
overload alias) versus a thin ``@overload`` + :keyword:`match` router
over shared primitives — is deliberately **parked on #261**, to be
settled *together* with the C / F / S core relocation
(:ref:`integral-kernel-category`). Once the iso / :math:`(n,2n)`
group-transfer, the :math:`R\circ\Lambda\circ M` reconstruction, and the
rank-1 production core move to shared ``transport/`` primitives, the
natural shape becomes a thin typed router over named single-sourced
primitives, and the ``match``-based form reads best (it also narrows away
the ``psi.bulk`` cast that #262 tracks). Deciding the spelling before the
cores move would be premature; the sharing should dictate the form. Until
then, Pattern M is correct, runtime-untouched, and pinned by the C6 gate
(``tests/sn/operators/test_operators_apply_typed.py`` — static
``assert_type`` pins with mutation-verified teeth, plus a runtime
dispatch-parity check on the aliased public ``apply``).


Pure-L streaming + the affine collision split
==============================================

The streaming leaf :math:`L` (:class:`~orpheus.sn.operators.streaming.StreamingOperator`)
computes **pure** :math:`\sigma`-free streaming directly: its ``apply``
is the named :math:`\sigma`-free
:meth:`~orpheus.sn.loss_representation.LossRepresentation.streaming_action`
leaf, the spatial streaming :math:`\Omega\cdot\nabla\psi` plus the
curvilinear angular redistribution, with NO collision diagonal.  The
collision diagonal :math:`C = M[\sigma_t]`
(a :class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`,
the §5.7 multiplier promotion — now literally :math:`C = M[\sigma_t]`) is
the separate shared leaf, and the composition
:math:`L + C` recovers the full within-group loss.

The discrete within-group WDD matvec is **affine in** :math:`\sigma` in
the forward direction:

.. (vv-status rationale) #257 S8b — the σ-free streaming primitive. The
   intrinsic σ-freedom of pure L (its apply reads no σ) is a SOFTWARE
   invariant pinned by the foundation catcher C1
   (``tests/sn/operators/test_pure_L_sigma_free.py``, with a Mode-11
   σ-leak mutation that reddens C1); the affine relation
   ``M(σ)ψ = streaming_action(ψ) + σ⊙ψ`` and the byte-identical (L+C)
   recovery are pinned by ``test_streaming_operator_decomposition.py``
   and ``test_loss_action_convention.py``. The streaming discretization
   is single-sourced through ``loss_action`` at σ = 0.
.. vv-status: streaming-action-pure-l documented

.. math::
   :label: streaming-action-pure-l

   M(\sigma)\,\psi \;=\; \underbrace{\text{streaming\_action}(\psi)}_{L\,\psi,
       \;\sigma\text{-free}} \;+\; \sigma_t \odot \psi
   \qquad\Longleftrightarrow\qquad
   \text{streaming\_action}(\psi) \;=\; \texttt{loss\_action}(0, \psi)

so :math:`L` reads no :math:`\sigma`: the curvilinear Carlson
coupled-pole seed's :math:`\sigma`-dependence is exactly the collision
diagonal it injects, which cancels into :math:`\sigma_t\odot\psi` and
belongs to :math:`C` (ERR-058 / #195 made the seed σ-independent, which
is what licenses the carve).  The streaming discretization lives ONCE in
``loss_action``; ``streaming_action`` is single-sourced from it at
:math:`\sigma = 0` (``coding-elegance`` Pattern 2), so there is no twin
σ-free walk.

Why the matvec is affine in :math:`\sigma`
------------------------------------------

The discrete within-group cell balance is the single source of the
affine structure. In the geometry-agnostic 1-D scan
(:meth:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.affine_scan_coefficients`)
and in the curvilinear cell update
(:func:`~orpheus.transport.spatial.diamond.cell_balance_for_streaming`) the
WDD cell-average solves

.. math::
   :label: streaming-action-cell-balance

   S\,\bar\psi \;=\; \underbrace{Q\,V \;+\; \text{(upstream-face inflow)}}_{\text{source} + \text{streaming numerator}},
   \qquad
   S \;=\; \underbrace{S_{\rm stream}}_{\text{geometric}}
       \;+\; \underbrace{\sigma_t\,V}_{\text{collision}},

where the cell-balance diagonal :math:`S` is the **sum** of a purely
geometric streaming term :math:`S_{\rm stream}` and the collision
volume term :math:`\sigma_t\,V`. In the production code
(:func:`~orpheus.transport.spatial.diamond.DiamondDifference._cartesian_streaming_diagonal`)
the Cartesian scan denominator is literally
``denom = reaction_xs + Σ_axes (2 g_axis)`` with ``reaction_xs`` the
collision term and ``2 g_axis = 2|μ_axis|/Δ_axis`` the geometric
streaming term; the curvilinear form
(:func:`~orpheus.transport.spatial.diamond.DiamondDifference.affine_scan_coefficients`)
is ``denom = geometric_streaming_term + collision_volume_term`` with
``geometric_streaming_term = 2|μ|·A_down + (ΔA/w)·c_out``. The
collision cross section enters the diagonal **purely additively** —
:math:`S` is affine in :math:`\sigma_t` with unit slope :math:`V`. That
is exactly :eq:`streaming-action-pure-l`: the forward matvec
:math:`M(\sigma)\psi` is the *application* of this diagonal (not its
inverse), so the collision contribution is the clean additive term
:math:`\sigma_t\odot\psi` that the pure-L leaf simply does not carry.

The subtle part is the **curvilinear angular closure**. For sphere /
cylinder the Carlson coupled-pole seed
(``precompute_psi_state``) builds a half-angle starting flux from the
flat moment :math:`\bar\phi_0` via a recursion whose coefficients —
:math:`\bar Q = \sigma_t\,\bar\phi_0` and a per-pole denominator
:math:`\Delta r\,\sigma_t + 2` — are themselves :math:`\sigma_t`-bearing.
At face value this looks like a :math:`\sigma_t`-dependence that lives
in the *streaming* (angular-redistribution) term, which would break the
clean additive split. It does not. The seed's :math:`\sigma_t`-dependence
is **exactly the collision diagonal it injects** into the half-angle
balance; when the cell update assembles
:math:`m_{\rm full} = (S\,\bar\psi - \text{numerator})/V`, that injected
collision exactly cancels the seed's :math:`\sigma_t`-bearing
contribution, leaving the redistribution term :math:`(1-\mu^2)/r\,
\partial\psi/\partial\mu` (the genuine angular streaming) **independent of
:math:`\sigma_t`**. The net :math:`\sigma_t` that survives is the single
:math:`\sigma_t\,V` collision term in :math:`S` — and nothing else.

ERR-058 / #195 is what *licenses* writing this down as a software
invariant. Before that fix, the curvilinear Carlson seed used a
:math:`\sigma_t`-coupled angular-edge extrapolation that left a residual
:math:`\sigma_t`-dependence in the redistribution; the decomposition
test's top docstring still carried the now-stale claim that
``matvec(σ=0)`` was "3–13 % wrong for curvilinear". ERR-058 replaced the
seed with a :math:`\sigma_t`-independent ``AngularEdgeExtrapolation``,
which made the curvilinear matvec genuinely affine in :math:`\sigma_t` —
and that is the precondition the pure-L carve depends on. (Issue #282
route (a) later retired that strategy *class*, but the
:math:`\sigma_t`-affinity invariant survives — the direct
starting-direction march and the inlined ``edge_extrapolated_seed`` are
both affine in :math:`\sigma_t` exactly like the bulk walk; see
:ref:`sn-282-direct-starting-direction-solve`.) The lesson
catalogued in :mod:`orpheus.sn.loss_representation` is to probe the live
behaviour rather than trust the prose: the affinity is an *empirical*
fact that must be re-verified, not a transcribed claim.

Probe evidence and the retired fold
------------------------------------

The carve does **not** duplicate the ~400-line discretization walk to
produce a separate :math:`\sigma`-free streaming kernel. ``loss_action``
is monolithic in :math:`\sigma` (the cross section is threaded into the
Cartesian ``residual_kernel_batch``, the curvilinear
``cell_balance_for_streaming``, *and* the Carlson seed
``precompute_psi_state``), so a hand-separated streaming walk would be a
twin path (Cardinal Rule 2 violation). Instead
:meth:`~orpheus.sn.loss_representation.LossRepresentation.streaming_action`
is **single-sourced** from the same walk at :math:`\sigma = 0`
(``coding-elegance`` Pattern 2 — name the primitive, do not clone the
algebra). Two in-process probes establish that this is value-correct,
not merely convenient:

.. list-table:: Pure-L σ-freedom — measured drift (#257 S8b)
   :header-rows: 1
   :widths: 46 18 18 18

   * - Probe
     - slab / Cartesian
     - sphere
     - cylinder
   * - Affine relation
       :math:`\texttt{loss\_action}(0,\psi) = \texttt{loss\_action}(\sigma,\psi)-\sigma\odot\psi`
       (bulk; boundary strict 0 ULP)
     - ≤ 32 ULP
     - ≤ 2 ULP
     - ≤ 72 ULP
   * - σ-leak test
       :math:`\texttt{streaming\_action}(\,\cdot\,;\,\sigma_a) = \texttt{streaming\_action}(\,\cdot\,;\,\sigma_b)`
       for two **wildly different** :math:`\sigma` fields
     - .. centered:: ≤ 64 ULP (rel ~1e-16), all geometries
     -
     -
   * - Pure-L leaf vs retired ``(L+C)−C`` fold
       (the matvec :math:`L\psi` the carve replaced)
     - ≤ 16 ULP
     - ≤ 16 ULP
     - ≤ 16 ULP

The σ-leak test is the decisive one: applying the streaming leaf to the
same flux under two completely different cross-section fields produces
byte-stable results (relative drift at the floating-point floor), so the
leaf demonstrably reads no :math:`\sigma`. The first row is the
quantitative form of the cancellation argument above — the difference
between the full loss and the pure stream is *exactly*
:math:`\sigma\odot\psi` to the floating-point floor, with the curvilinear
geometries (where the cancellation runs through the angular closure) the
tightest of all (sphere ≤ 2 ULP). The third row records that the carve
re-associates the floating-point reduction tree relative to the retired
``(L+C)−C`` fold (the old way of obtaining :math:`L\psi`: build the full
loss, subtract the collision diagonal) — the values agree to ≤ 16 ULP,
well inside the dimensionally-explainable single-step bound (per
``vv-principles`` § "Bit-identity vs principled-equivalence", criterion 3).

The carve is therefore **principled-equivalent, not bit-identical**, on
the *streaming-leaf* matvec — and **byte-identical** on the
:math:`(L+C)` composite matvec and the WDD sweep, which were not touched
(:meth:`~orpheus.sn.operators.streaming.InvertibleOperator.apply` still computes
:math:`M(\sigma_t)\psi` through the same ``loss_action`` call). The
software invariant "pure :math:`L` reads no :math:`\sigma`" is pinned by
the foundation catcher
:func:`tests.sn.operators.test_pure_L_sigma_free.test_c1_pure_L_apply_is_sigma_free`,
which carries a Mode-11 σ-leak mutation (monkeypatch a σ-leaking
``streaming_action`` stub) that reddens the gate — so the gate is
verified to be *able* to see a regression, not merely green. The affine
relation and the byte-identical :math:`(L+C)` recovery are pinned by
``test_streaming_operator_decomposition.py`` and
``test_loss_action_convention.py``.

This affine structure is the algebraic foundation of the next subsection:
because the *forward* application is affine in :math:`\sigma` (additive,
distributive), the leaves compose additively — but the *inverse* is not,
which is precisely why ``solve`` cannot live on the leaves.


.. _apply-solve-asymmetry:

apply is linear in the operator; solve is not
=============================================

The single most important structural fact about the
:math:`L+C` algebra is an **asymmetry between the two primitive actions**
of :ref:`Definitions <operator-algebra>`: forward application
(:eq:`operator-apply`) is *linear in the operator*, but inversion
(:eq:`operator-solve`) is *not*. This is why ``apply`` and ``solve`` are
**two faithful views of the same operator only for the bundled**
:class:`~orpheus.sn.operators.streaming.InvertibleOperator` (:math:`= L+C`), and
**never** for the individual streaming / collision leaves. It is the
mathematical content behind the :ref:`three-layer operator surface
<capability-set-semantics>`: ``apply`` lives on the leaves; ``solve``
lives on the bundle.

The asymmetry
-------------

Forward application **distributes over the operator sum**, because
applying a sum of linear maps is the sum of the applications:

.. math::
   :label: apply-distributes

   (L + C)\,\psi \;=\; L\,\psi \;+\; C\,\psi.

Inversion does **not** distribute:

.. math::
   :label: solve-does-not-distribute

   (L + C)^{-1} \;\neq\; L^{-1} \;+\; C^{-1}.

The inverse is a :math:`1/x`-shaped functional of the operator, and
:math:`1/x` does not distribute over :math:`+`. The crispest sanity
anchor is the scalar case: :math:`(3+5)^{-1} = 1/8 = 0.125`, whereas
:math:`3^{-1} + 5^{-1} = 0.533`. The two are not equal, not even close.
``L.solve(q) + C.solve(q)`` would be the operator version of
:math:`q/L + q/C`, which is emphatically **not** :math:`q/(L+C)`. This is
the same fact that the :ref:`composition algebra
<composition-algebra>` table encodes as "``OperatorSum`` does not
propagate ``solve``": no general algorithm exists for :math:`(A+B)^{-1}`
from :math:`A^{-1}` and :math:`B^{-1}` (Sherman–Morrison–Woodbury applies
only under low-rank structure, which the dense collision diagonal is not).

What each separate inverse would mean physically
------------------------------------------------

The within-group transport balance is the operator equation
:math:`(L+C)\,\psi = q`, i.e.

.. math::
   :label: apply-solve-within-group-balance

   \Omega\cdot\nabla\psi \;+\; \Sigt{}\,\psi \;=\; q,

with :math:`L = \Omega\cdot\nabla` the streaming (advection) operator and
:math:`C = M[\Sigt{}]` the collision diagonal. Each *separate* inverse
solves a *different, decoupled* problem:

.. list-table:: The three inverses solve three different problems
   :header-rows: 1
   :widths: 16 30 54

   * - Inverse
     - Equation it solves
     - Physical meaning
   * - :math:`L^{-1}`
     - :math:`\Omega\cdot\nabla\psi = q`
     - **Pure advection, no absorption** — the flux if neutrons streamed
       freely and never collided. A formal inverse only: pure streaming
       is rank-deficient (see below).
   * - :math:`C^{-1}`
     - :math:`\Sigt{}\,\psi = q \;\Rightarrow\; \psi = q/\Sigt{}`
     - **Infinite-medium / no-leakage flux** — purely local, the flux at
       a point set entirely by the local collision rate, with no spatial
       coupling whatsoever.
   * - :math:`(L+C)^{-1}`
     - :math:`\Omega\cdot\nabla\psi + \Sigt{}\,\psi = q`
     - **The coupled balance** — the flux that satisfies *both* loss
       mechanisms simultaneously, everywhere. This is the true transport
       solution, and it is **not** the sum of the two decoupled problems.

The true solution :math:`(L+C)^{-1}q` is genuinely coupled: streaming
moves a neutron from cell to cell while collision removes it, and the
two compete at every point. Solving them separately and adding throws
away the competition. That is *why* a Neumann-style series is the only
honest way to express the coupled inverse through the parts (below).

.. note::

   ``L.solve`` is **not a live call** in ORPHEUS. The streaming leaf
   :class:`~orpheus.sn.operators.streaming.StreamingOperator` is
   adjointable but **not invertible** (it exposes ``apply`` and
   ``apply_transpose``, reports
   :attr:`~orpheus.numerics.operator.LinearOperator.is_invertible`
   ``= False``, and declares **no** ``inverse()`` / ``solve`` — the
   *structural* arm of the two-kinds split, Design C) — precisely
   because pure streaming is rank-deficient (without a collision term
   the within-group cell balance is singular;
   :eq:`streaming-action-cell-balance` has :math:`\sigma_t\,V = 0`, so
   :math:`S` degenerates to the geometric streaming term alone, which has
   a zero-mode the inflow boundary condition cannot pin in general). The
   :math:`L^{-1}` row above is therefore the **mathematical** advection
   inverse — "even if :math:`L` alone were inverted, it would solve pure
   advection" — not a method you can invoke (asking for it is a *static*
   error, not a runtime raise). The collision leaf
   :math:`C = M[\sigma_t]`
   (a :class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`)
   *does* report
   :attr:`~orpheus.numerics.operator.LinearOperator.is_invertible`
   ``= True`` and carry ``solve``, but **only** when
   :math:`\min|\sigma| > 0` (the multiplier spectrum law
   :math:`\mathrm{spec}(M[\sigma]) = \mathrm{ess\,range}(\sigma)`) — it
   is the *value-dependent* arm: it always declares ``inverse()`` /
   ``solve`` and refuses eagerly with
   :class:`~orpheus.numerics.operator.NotInvertible` on a zero
   coefficient. When invertible, ``C.solve(q) = q/\sigma``
   is the infinite-medium flux of the :math:`C^{-1}` row, computed as an
   element-wise division.

The crispest proof — the WDD cell denominator
---------------------------------------------

The non-separability is visible inside a **single cell**, before any
global coupling enters. The diamond-difference cell update solves the
balance :eq:`streaming-action-cell-balance` for the cell-average flux:

.. math::
   :label: apply-solve-cell-resolvent

   \bar\psi
   \;=\;
   \frac{Q\,V \;+\; \text{(upstream-face inflow)}}{S},
   \qquad
   S \;=\; S_{\rm stream} \;+\; \sigma_t\,V,

so the resolvent (the per-cell ``solve``) **divides by the SUM**
:math:`S = S_{\rm stream} + \sigma_t\,V`. In the production code this is
literally ``inverse_denom = 1.0 / denom`` with ``denom = streaming_term
+ collision_term``
(:func:`~orpheus.transport.spatial.diamond.DiamondDifference.affine_scan_coefficients`,
:func:`~orpheus.transport.spatial.diamond.DiamondDifference.cartesian_scan_coefficients`).
Now compare the three inverses on this single cell:

- :math:`L^{-1}` would divide by :math:`S_{\rm stream}` **alone**
  (set :math:`\sigma_t = 0`) — the pure-streaming denominator.
- :math:`C^{-1}` would divide by :math:`\sigma_t\,V` **alone** (no
  upstream-face coupling at all) — the local denominator.
- The coupled resolvent divides by :math:`S_{\rm stream} + \sigma_t\,V`,
  and

  .. math::
     :label: apply-solve-denominator-inequality

     \frac{1}{S_{\rm stream} + \sigma_t\,V}
     \;\neq\;
     \frac{1}{S_{\rm stream}} \;+\; \frac{1}{\sigma_t\,V}.

The two losses are **added before the division** (the additive structure
of :eq:`streaming-action-pure-l`); you must invert the *sum*. This is the
single-cell shadow of the global inequality
:eq:`solve-does-not-distribute`. The forward matvec, by contrast,
*multiplies* by the cell-balance diagonal :math:`S` — and multiplication
by a sum distributes (:eq:`apply-distributes`), which is exactly why
``apply`` survives on the leaves while ``solve`` does not.

The Neumann series — the only honest way through the parts
----------------------------------------------------------

If one insists on expressing :math:`(L+C)^{-1}` through the individual
inverses, the only correct expression is an **infinite operator-splitting
(Neumann) series**, not a finite sum. Splitting around the collision
diagonal :math:`C` (which is the cheap-to-invert leaf, ``C.solve(q) =
q/\sigma``):

.. math::
   :label: apply-solve-neumann-series

   (L+C)^{-1}
   \;=\;
   \bigl[C\,(I + C^{-1}L)\bigr]^{-1}
   \;=\;
   (I + C^{-1}L)^{-1}\,C^{-1}
   \;=\;
   \sum_{k=0}^{\infty} (-1)^{k}\,(C^{-1}L)^{k}\,C^{-1},

i.e.

.. math::
   :label: apply-solve-neumann-expansion

   (L+C)^{-1}
   \;=\;
   C^{-1} \;-\; C^{-1}L\,C^{-1} \;+\; C^{-1}L\,C^{-1}L\,C^{-1} \;-\;\cdots,

which converges when the spectral radius
:math:`\rho(C^{-1}L) < 1`. The leading term :math:`C^{-1}` is the
infinite-medium flux; every subsequent term is a streaming correction.
A *finite* truncation — and in particular the one-term sum
:math:`L^{-1} + C^{-1}` — is **never** the coupled inverse. The closest
clean closed form involving both inverses is the **parallel** (resistors
-in-parallel) identity

.. math::
   :label: apply-solve-parallel-identity

   \bigl(L^{-1} + C^{-1}\bigr)^{-1}
   \;=\;
   L\,(L+C)^{-1}\,C,

which is still **not** :math:`(L+C)^{-1}` — it is the harmonic
combination, related to but distinct from the coupled inverse.

This is more than an algebraic curiosity: the transport-native Neumann
series is the **source-iteration / collision-number expansion** itself.
The full within-group-plus-scattering problem
:math:`(L+C-S)\psi = q` is solved by splitting off the scattering source
:math:`S` and summing the series

.. math::
   :label: apply-solve-source-iteration-series

   \psi
   \;=\;
   \sum_{k=0}^{\infty}
     \bigl[(L+C)^{-1}S\bigr]^{k}\,(L+C)^{-1}\,q,

where the **sweep** :math:`(L+C)^{-1}` is the per-term inverter and the
outer source iteration sums the series. This is the Peierls
collision-number expansion (each term :math:`k` is the flux of neutrons
that have scattered exactly :math:`k` times). The series converges with
spectral radius :math:`\rho\bigl[(L+C)^{-1}S\bigr] \le \max_g
\Sigma_{s,g}/\Sigma_{t,g} = c` (the scattering ratio;
:ref:`affine-typed-field-algebra` documents the matching contraction
:math:`M = (L+C)^{-1}(S+B)` carried as a typed diagnostic on the SI
iterate). The sweep :math:`(L+C)^{-1}` being a *single bundled* inverse —
not :math:`L^{-1} + C^{-1}` — is exactly the point: it is the WDD
forward-substitution on :eq:`apply-solve-cell-resolvent`, dividing by the
summed denominator cell-by-cell in inflow-to-outflow order. See Lewis &
Miller, *Computational Methods of Neutron Transport* ([LewisMiller1984]_,
§3.2 for the sweep as the discrete-ordinates resolvent and §4 for the
source-iteration / Neumann scattering series), and Adams & Larsen 2002
([AdamsLarsen2002]_, §II for the spectral radius :math:`\rho = c`).

Why this is the right architecture, not a limitation
----------------------------------------------------

Invertibility is a property of the **sum**, not of the parts. That is
exactly why :math:`L + C` is packaged as one
:class:`~orpheus.sn.operators.streaming.InvertibleOperator`: the
:class:`~orpheus.numerics.operator.OperatorSum` that *carries* the
WDD sweep as its ``.solve``. The asymmetry maps cleanly onto the two
sides of the algebra:

- **apply lives on the faithful separate leaves.** Pure streaming
  :math:`L` (:class:`~orpheus.sn.operators.streaming.StreamingOperator`, the
  :math:`\sigma`-free :eq:`streaming-action-pure-l` leaf) and collision
  :math:`C = M[\sigma_t]`
  (a :class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`)
  each advertise
  ``apply``, and their applications compose additively
  (:eq:`apply-distributes`). The forward direction is affine in
  :math:`\sigma` (the previous subsection), so the leaf decomposition is
  *faithful*: :math:`(L+C)\psi = L\psi + C\psi` holds exactly.
- **solve belongs to the bundled unit.** Only
  :class:`~orpheus.sn.operators.streaming.InvertibleOperator` reports
  :attr:`~orpheus.numerics.operator.LinearOperator.is_invertible`
  ``= True`` and carries a direct-sweep ``solve``; the leaves do not
  (streaming has no ``solve`` at all; collision's ``solve`` is the
  *local* :math:`q/\sigma`, which is the
  :math:`C^{-1}` of a *different* problem, never the coupled inverse).
  The :class:`~orpheus.numerics.operator.OperatorSum` deliberately
  **does not propagate** ``solve`` (:ref:`composition-algebra`); the
  :class:`~orpheus.sn.operators.streaming.InvertibleOperator` *adds it back* via the
  SN-specific algebraic identity "WDD sweep :math:`\approx (L+C)^{-1}`"
  ([LewisMiller1984]_ §3.2). The composite owns ``apply``, ``solve``, and
  ``apply_transpose`` as three actions of **one** operator on a single
  shared
  :class:`~orpheus.sn.loss_representation.LossRepresentation` (L21 —
  "matvec ≡ sweep"); :meth:`InvertibleOperator.apply` and
  :meth:`InvertibleOperator.solve` single-source :math:`\sigma` from the
  collision diagonal, so they cannot disagree on which loss they invert.

The :ref:`three-layer operator surface <capability-set-semantics>` is
what makes this architecture *enforced* rather than merely *intended*: a
downstream Krylov consumer that asks for the inverse of a bare
:class:`~orpheus.sn.operators.streaming.StreamingOperator` cannot even
spell ``L.inverse()`` (the streaming leaf declares no such method — a
*static* error), and a generic sum that has not been promoted to
:class:`~orpheus.sn.operators.streaming.InvertibleOperator` returns the
iterative splitting :class:`~orpheus.numerics.green_operator.GreenOperator`
rather than the direct sweep — never silently handed
:math:`L^{-1} + C^{-1}` (a meaningless answer to a problem nobody
posed). The asymmetry between :eq:`apply-distributes` and
:eq:`solve-does-not-distribute` is, in this sense, the *reason
invertibility is a per-instance predicate on the operator, not a flag in
a parallel string registry*.

.. vv-status: apply-distributes documented
.. vv-status: solve-does-not-distribute documented


.. _capability-set-semantics:

The three-layer operator surface
================================

Many transport operators have **no** efficient inverse action, and this
is a load-bearing fact, not an inconvenience:

- The scattering source operator :math:`S` has rank in the thousands
  of unknowns and is never inverted directly (the source iteration
  scheme exists precisely to avoid that inversion).
- The fission source operator :math:`F` is rank-deficient — it
  projects onto the fission spectrum :math:`\chi`. There is no
  inverse.
- A Jacobi-preconditioned matvec wrapped for BiCGSTAB has ``apply``
  but no ``solve`` (the solve is the iterative scheme itself, not
  this operator).

Subclassing or abstract methods would force these classes to provide
``solve`` stubs that raise :class:`NotImplementedError`. That is the
**harmful-stub anti-pattern**: downstream Krylov consumers would only
discover the failure mid-iteration. The honest surface for a missing
ability is **method absence**, not an advertising flag — but *how* that
honesty is enforced is the subject of this section.

.. note:: **The retired capability frozenset (#226 taxonomy step 6).**

   Through step 5 the advertisement was a stringly-typed
   ``capabilities: frozenset[str]`` class property listing the ``CAP_APPLY``
   / ``CAP_SOLVE`` / ``CAP_APPLY_TRANSPOSE`` tags an operator supported,
   plus a ``MissingCapability`` exception raised when a composition asked
   for an absent tag. **Carve P4 retired it from every operator** — leaves,
   composers, aggregators, and shims. The frozenset was a **parallel
   registry that could silently drift from the actual method surface**: a
   class could advertise ``CAP_SOLVE`` yet have a broken ``solve`` (or, as
   the collision leaf once did, a ``solve`` that produced a silent NaN), and
   nothing forced the string set to track the code. The replacement makes
   the **surface itself the single source of truth**: an ability exists
   exactly where its method exists, and a runtime predicate reads the live
   structure and values rather than a cached string. The design below is
   the user-locked "Design C + B" (2026-07-02).

The three layers, per axis
--------------------------

Each optional axis — **inverse** and **adjoint** — now has exactly three
layers, and each layer carries the truth it *alone* can express:

.. list-table::
   :header-rows: 1
   :widths: 20 34 46

   * - Layer
     - Inverse axis
     - Adjoint axis
   * - **Predicate** (runtime truth)
     - :attr:`~orpheus.numerics.operator.LinearOperator.is_invertible`
     - :attr:`~orpheus.numerics.operator.LinearOperator.is_adjointable`
   * - **Operator-returning method**
     - ``inverse()`` — declared *per-class*
     - :attr:`~orpheus.numerics.operator.LinearOperator.H` — hosted on the base
   * - **Realization verb**
     - ``solve``
     - ``apply_transpose``

- The **predicate** is the runtime, instance-accurate truth. Unlike an
  ``isinstance`` check — which sees only class-level method presence — it
  reads the operator's actual structure AND values: a zero-coefficient
  :class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`
  reports ``is_invertible = False``; a sum reports its *leading* term
  (the left-spine head — :ref:`green-operator`); a composite derives its
  value recursively from the operands, never a cached registry. The
  default is ``False`` — an operator is invertible or adjointable only by
  explicit override.
- The **operator-returning method** is the canonical act: it returns an
  *operator*, not a vector. ``inverse()`` returns the inverse operator
  (a member of the inverse family — the sweep, the
  :class:`~orpheus.numerics.green_operator.GreenOperator`, the
  :class:`~orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator`);
  :attr:`~orpheus.numerics.operator.LinearOperator.H` returns the
  metric-correct Hilbert adjoint wrapper.
- The **realization verb** is the raw numerical act — ``solve`` maps a
  right-hand side to a solution vector, ``apply_transpose`` applies the
  Euclidean transpose — present **exactly where a native realization
  exists**, never as an exists-but-raises stub. The wrap-delegate
  inverse family delegates through ``solve``; the composer transpose
  laws recurse through ``apply_transpose``.

.. _design-c-structural-value-split:

Design C — the structural-vs-value split
----------------------------------------

The load-bearing insight is that **two mathematically distinct kinds of
non-invertibility deserve two honest surfaces**, and conflating them was
a false dichotomy (see :ref:`design-c-false-dichotomy` below). The two
kinds:

.. list-table::
   :header-rows: 1
   :widths: 22 30 48

   * - Kind
     - The claim
     - Surface
   * - **Structural**
     - *This TYPE has no inverse* — the map is mathematically
       non-invertible for every instance.
     - ``inverse()`` is **not declared at all**. Asking for it is a
       *static* error (pyright reports ``reportAttributeAccessIssue`` at
       the call site); at runtime the attribute is simply absent.
   * - **Value-dependent**
     - *This TYPE supports inversion, this INSTANCE refuses* — a
       zero-coefficient diagonal, a sum with a non-invertible leading
       term, a product with a singular factor.
     - ``inverse()`` **is declared** and raises
       :class:`~orpheus.numerics.operator.NotInvertible` **eagerly**, at
       construction of the inverse and never mid-iteration.

The structural leaves are
:class:`~orpheus.numerics.operator.ZeroOperator`, the incoming-ordinate
and periodic masks, the source dyads
(:class:`~orpheus.numerics.operator.RankOneOperator`), and the
transport source leaves —
:class:`~orpheus.sn.operators.streaming.StreamingOperator` (:math:`L`),
:class:`~orpheus.transport.operators.scattering.ScatteringOperator`
(:math:`S`),
:class:`~orpheus.transport.operators.fission.FissionOperator`
(:math:`F`), and the boundary operator
:class:`~orpheus.sn.operators.boundary.SNBoundaryOperator` (:math:`B`).
For each, ``op.inverse()`` does not type-check — the absence *is* the
honest surface, and the compiler enforces it.

The value-dependent operators are
:class:`~orpheus.numerics.operator.DiagonalOperator`,
:class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`,
and the composers when their invertibility fails a value test
(:class:`~orpheus.numerics.operator.OperatorSum` with a non-invertible
head, :class:`~orpheus.numerics.operator.OperatorProduct` with a
singular factor,
:class:`~orpheus.numerics.operator.ScaledOperator` of a non-invertible
operand). They *always* declare ``inverse()`` (the type *can* invert),
and refuse loudly and eagerly when the specific values do not permit it.

.. _typeguard-bridge:

The checked bridge — a PEP-647 ``TypeGuard``
--------------------------------------------

The predicate is a *runtime* fact; ``inverse()`` living only on some
classes is a *static* fact. The bridge that converts one into the other
is a pair of free functions,
:func:`~orpheus.numerics.operator.invertible` and
:func:`~orpheus.numerics.operator.adjointable`, each a PEP-647
``TypeGuard`` narrowing a
:class:`~orpheus.numerics.operator.LinearOperator` to the
:class:`~orpheus.numerics.operator.SupportsInverse` /
:class:`~orpheus.numerics.operator.SupportsAdjoint` narrowing target:

.. code-block:: python

   if invertible(op):          # runtime check on op.is_invertible
       y = op.inverse().apply(b)   # pyright now permits .inverse() — no cast

**The runtime check IS the static permission.** You cannot obtain the
permission without executing the check, and deleting a guard un-narrows
the call so CLI pyright REDs — the guards are *type-load-bearing* (this
is verified as a static tooth: spec §39.1). This finally fulfils the
original charter of the
:class:`~orpheus.numerics.operator.SupportsInverse` /
:class:`~orpheus.numerics.operator.SupportsAdjoint` Protocols (a static
contract backed by a runtime property) through a *checked* bridge — the
carve deleted all four ``cast(SupportsInverse, …)`` sites and all ten
``solve`` / ``apply_transpose`` ``# type: ignore`` comments they used to
require.

Two subtleties fix the exact construct:

- **``TypeGuard``, deliberately NOT ``TypeIs``.** The predicate is
  value-dependent: a zero-coefficient multiplier structurally *has*
  ``inverse()`` while reporting ``is_invertible = False``. Only the
  *one-directional* promise is honest — ``True`` licenses the call;
  ``False`` makes no static claim (a ``TypeIs`` would wrongly narrow the
  negative branch too, asserting the operand is *not* a
  ``SupportsInverse`` when structurally it may still be one).
- **A free function, not a method.** PEP 647 narrowing applies only
  through a call expression, and a method form narrows its first
  *explicit* argument, never ``self`` — so there is no
  ``op.is_invertible``-style property spelling that could narrow the
  operand. The narrowing target
  :class:`~orpheus.numerics.operator.SupportsInverse` was *promoted* to
  **extend** :class:`~orpheus.numerics.operator.LinearOperator`, so a
  guarded branch keeps the whole algebra (``apply``,
  :attr:`~orpheus.numerics.operator.LinearOperator.H`, the composition
  dunders) alongside the licensed ``inverse()``.

.. _base-hosting-rule:

The base-hosting rule
---------------------

**A method lives on the base** :class:`~orpheus.numerics.operator.LinearOperator`
**iff a universal realization exists.**

- :attr:`~orpheus.numerics.operator.LinearOperator.H` **is base-hosted**:
  the Hilbert adjoint has one generic realization — the
  ``_AdjointOperator`` wrapper that applies the metric once and delegates
  to ``apply_transpose`` — valid for *any* adjointable operator, so
  ``.H`` is defined once on the base with an eager
  :class:`~orpheus.numerics.operator.MissingAdjoint` gate.
- ``inverse()``, ``solve``, and ``apply_transpose`` are **NOT
  base-hosted**: there is no universal inverse (each structure inverts
  differently — a diagonal by reciprocal, a triangular sweep by
  forward-substitution, a full loss by a Neumann splitting), no universal
  transpose realization, and no universal solve. Each lives per-class,
  exactly where its realization exists.

This is why the structural non-invertibility surface *works*: because
``inverse()`` is not on the base, a class that omits it genuinely has no
such attribute, and ``ZeroOperator().inverse()`` is a static error
rather than a stub that raises.

.. _design-b-native-solve:

Design B — ``solve`` pruned to native realizations
--------------------------------------------------

The realization verb ``solve`` is now present **only where a native
realization exists** — never as an exists-but-raises stub, and never
duplicating what ``.inverse().apply`` already does. This executes the
"one public surface = predicate + operator-returning method" ruling
(taxonomy §11): *solving with an operator IS applying its inverse
object*, ``A.inverse().apply(b)``.

**Deleted** (the algebra-closed and driver-realized kinds):

- :class:`~orpheus.numerics.operator.OperatorSum` — a generic sum's
  inverse action is driver-realized by the
  :class:`~orpheus.numerics.green_operator.GreenOperator` (the
  preconditioned splitting), not a substrate verb. (The
  sweep-invertible ``(L+C)``
  :class:`~orpheus.sn.operators.streaming.InvertibleOperator` subclass
  keeps its own direct-sweep ``solve`` — that IS a native realization.)
- :class:`~orpheus.numerics.operator.IdentityOperator`,
  :class:`~orpheus.numerics.operator.PermutationOperator`,
  :class:`~orpheus.numerics.operator.ScaledOperator`,
  :class:`~orpheus.numerics.operator.TensorProductOperator` — the
  **algebra-closed** kinds, whose inverse is itself a first-class
  *forward* operator (a permutation's inverse is a permutation, a
  scaling's is a scaling): solving is just ``.inverse().apply``.
- the reflective boundary shim's forward.

**Kept** (the native-realization face, what the wrap-delegate inverse
family wraps): :class:`~orpheus.numerics.operator.DiagonalOperator`
(with its value guard, now raising
:class:`~orpheus.numerics.operator.NotInvertible`),
:class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`,
the sweep composites, the mixin's un-invert, and
:class:`~orpheus.numerics.operator.OperatorProduct` — whose ``solve``
**re-routes** through each factor's canonical surface:

.. (vv-status rationale) Representational identity — the factor-wise
   ``OperatorProduct.solve`` re-route (Design B), NOT a solver claim
   (no eigenvalue / flux). The verifiable content is the §40 re-route
   gate battery in ``tests/numerics/test_operator.py`` (the dense
   ``np.linalg.solve`` anchor + the five-row factor-kind matrix vs the
   pre-carve baseline + the Mode-11 execution sentinel) and the §13 I2
   functoriality gate ``(AB)^{-1} = B^{-1}A^{-1}``.
.. vv-status: product-solve-reroute documented

.. math::
   :label: product-solve-reroute

   (A\,B)^{-1} b
   \;=\;
   B^{-1}\bigl(A^{-1} b\bigr)
   \;=\;
   \texttt{self.b.inverse().apply(self.a.inverse().apply(b))} .

The re-route is **bit-identical per factor kind** — each inverse object
delegates to the same realization the factor's own ``solve`` used to —
AND **total over the solve-retired kinds**: a permutation / scaling /
Green-invertible-sum factor, whose own ``solve`` Design B deleted, now
composes through its ``.inverse().apply`` (a permutation's inverse is a
first-class forward; a sum-in-a-product works via its
:class:`~orpheus.numerics.green_operator.GreenOperator`). The re-route is
gated at spec §40 by a structurally-independent dense
``np.linalg.solve(as_matrix, q)`` anchor plus a five-row factor-kind
matrix against a pre-carve baseline snapshot — with a Mode-11
counter-sentinel proving ``OperatorProduct.solve`` actually *executes*
under ``(A@B).inverse().apply`` (a value gate spelled ``b.inverse().apply
(a.inverse().apply(q))`` on both sides would be tautological).

The exception successors
------------------------

``MissingCapability`` split into **two** ``TypeError`` successors, one
per axis:

- :class:`~orpheus.numerics.operator.NotInvertible` — the inverse-axis
  refusal, raised eagerly by ``inverse()`` overrides (the
  value-dependent arm).
- :class:`~orpheus.numerics.operator.MissingAdjoint` — the adjoint-axis
  refusal, raised eagerly by
  :meth:`~orpheus.numerics.operator.LinearOperator.adjoint` /
  :attr:`~orpheus.numerics.operator.LinearOperator.H` and by the
  composer ``apply_transpose`` law bodies.

Both parent to :class:`TypeError`, carrying the retired
``MissingCapability``'s public contract forward: **no ``except`` clause
written against the old gate changes meaning**. (The migration was
staged: at wave W1 ``NotInvertible`` was born as a ``MissingCapability``
subclass so every landed ``pytest.raises(MissingCapability)`` stayed
green by inheritance while the new keystone went live; wave W2
re-parented both to ``TypeError`` and deleted the old class.)

.. _eager-adjoint-behavior-change:

The one behavior change — eager ``.H``
--------------------------------------

The carve is otherwise behavior-preserving; it has **exactly one**
observable behavior change (spec §38). Before, ``A.H`` on a
non-adjointable ``A`` *succeeded* — it constructed the wrapper
unconditionally — and the refusal was **lazy**, deferred to the
wrapper's first ``.apply``. Now
:meth:`~orpheus.numerics.operator.LinearOperator.adjoint` raises
:class:`~orpheus.numerics.operator.MissingAdjoint` **eagerly at
construction**:

.. code-block:: python

   def adjoint(self):
       if not adjointable(self):
           raise MissingAdjoint(...)   # eager — was lazy at .apply
       return _AdjointOperator(self)

A wrapper that could only fail at its first ``.apply`` is precisely the
broken-stub anti-pattern this module refuses; the
:func:`~orpheus.numerics.operator.adjointable` guard doubles as the
static bridge (the wrapper's constructor consumes the narrowed
:class:`~orpheus.numerics.operator.SupportsAdjoint`).

.. _design-c-false-dichotomy:

What was tried and rejected — the per-class-casts vs base-declaration false dichotomy
-------------------------------------------------------------------------------------

The design that Design C *replaced* framed the choice as a dichotomy
between two unappealing options, and the resolution was to recognise the
dichotomy as false (taxonomy §16 record):

- **Option A — declare ``inverse()`` on the base.** Then every call site
  type-checks, but a
  :class:`~orpheus.numerics.operator.ZeroOperator` (which mathematically
  has no inverse) would inherit a method it cannot honour — demoting the
  compiler's ability to catch ``Zero.inverse()`` misuse from a *static
  error* to a *runtime raise*. The honest "this type has no inverse"
  signal is lost.
- **Option B — keep ``inverse()`` per-class and ``cast`` at every call
  site.** Then structural absence is preserved (a static error on
  ``Zero.inverse()``), but every composer body that calls ``op.inverse()``
  on a ``LinearOperator``-typed operand needs a ``cast(SupportsInverse,
  op)`` and a ``# type: ignore`` — an *unchecked* assertion that could
  drift from the runtime truth exactly as the frozenset did.

The false premise was that *structural* absence and *value-dependent*
refusal are the same phenomenon. They are not: the
:class:`~orpheus.numerics.operator.ZeroOperator` case is structural
(Option A's failure), the zero-coefficient
:class:`~orpheus.numerics.operator.DiagonalOperator` case is
value-dependent (a method that *should* exist and refuse). **Design C
splits them** — structural leaves omit the method (Option B's win,
static error), value-dependent operators declare it and raise
``NotInvertible`` (Option A's ergonomics, honest runtime refusal) — and
replaces the *unchecked* ``cast`` (Option B's cost) with the *checked*
``TypeGuard`` bridge, which certifies the narrowing at runtime. The two
false horns dissolve because the two phenomena were never one.

.. note:: **What the static layer can and cannot certify.**

   Do NOT annotate a parameter with
   :class:`~orpheus.numerics.operator.SupportsInverse` to *demand*
   invertibility: the static layer can certify only **spelling** (the
   method exists on the class), never **solvability** (the value-level
   predicate). A zero-coefficient
   :class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`
   satisfies ``SupportsInverse`` structurally yet is not invertible.
   Guard with :func:`~orpheus.numerics.operator.invertible` at the
   ``LinearOperator``-typed call site instead — and only there, since a
   ``TypeGuard`` **replaces** (does not intersect) the declared type, so
   guarding an already-concrete operand would widen it.


.. _composition-algebra:

Composition algebra
===================

The composers in :mod:`orpheus.numerics.operator` implement the
following closure laws:

.. list-table::
   :header-rows: 1
   :widths: 22 26 26 26

   * - Composer
     - ``apply``
     - ``solve``
     - ``apply_transpose``
   * - :class:`~orpheus.numerics.operator.OperatorSum`
       (:math:`A + B`)
     - Both must have ``apply``.
     - **Does not propagate.** No general algorithm exists for
       :math:`(A + B)^{-1}` from :math:`A^{-1}, B^{-1}`
       (Sherman-Morrison-Woodbury applies only under low-rank
       structure).
     - Both must have ``apply_transpose``;
       :math:`(A + B)^T = A^T + B^T`.
   * - :class:`~orpheus.numerics.operator.OperatorProduct`
       (:math:`A\,B`)
     - Both must have ``apply`` (function composition).
     - Both must have ``solve``, **applied in REVERSE order**:
       :math:`(A\,B)^{-1} = B^{-1}\,A^{-1}`.
     - Both must have ``apply_transpose``, applied in REVERSE order:
       :math:`(A\,B)^T = B^T\,A^T`.
   * - :class:`~orpheus.numerics.operator.ScaledOperator`
       (:math:`\alpha\,L`)
     - Always preserved.
     - Preserved with division:
       :math:`(\alpha L)^{-1} = (1/\alpha)\,L^{-1}`. Zero scalar is
       rejected — use :class:`ZeroOperator` instead.
     - Preserved (scalars commute with transpose).
   * - :class:`~orpheus.numerics.operator.IdentityOperator`
       (:math:`I`)
     - Trivially yes.
     - ``solve`` is the same code path as ``apply``:
       :math:`I^{-1} = I`.
     - Trivially yes: :math:`I^T = I`.
   * - :class:`~orpheus.numerics.operator.ZeroOperator`
       (:math:`0`)
     - Trivially yes (returns ``np.zeros_like(x)``).
     - **Does not propagate.** The zero operator is not invertible.
     - Trivially yes (returns zero).

The closure rules above are advertised by the recursive
:attr:`~orpheus.numerics.operator.LinearOperator.is_invertible` /
:attr:`~orpheus.numerics.operator.LinearOperator.is_adjointable`
predicates and enforced by **eager guards** at composition time (an
``apply`` mismatch raises ``TypeError``, a non-adjointable ``.H`` raises
:class:`~orpheus.numerics.operator.MissingAdjoint`, a value-dependent
``inverse()`` raises :class:`~orpheus.numerics.operator.NotInvertible`),
never mid-iteration. A generic
:class:`~orpheus.numerics.operator.OperatorSum` carries no ``solve``
verb at all (Design B): solving with it IS applying its inverse object,
``.inverse().apply`` (the
:class:`~orpheus.numerics.green_operator.GreenOperator` splitting), so a
downstream consumer that spells :math:`L = A + 0` and asks for its
inverse action cannot silently receive :math:`A^{-1} + 0^{-1}`.


Cross-solver consumption (forward reference)
============================================

The Phase 0 algebra is the foundation for the cross-solver migration
sequence: SN, MoC, CP, and diffusion will each expose their natural
operator triple :math:`(A, S, F)` as
:class:`~orpheus.numerics.operator.LinearOperator` instances, and a
unified iteration driver in :mod:`orpheus.numerics.iteration` will
consume those triples without knowing which transport discretisation
they came from. That driver is **to be installed** under the Wave H
refactor (Issue 14). The Phase 0 module ships the algebra; subsequent
phases mount the producers and the consumer onto it.

Today the existing
:func:`build_transport_linear_operator <orpheus.sn.operators.streaming.build_transport_linear_operator>`
functions in :mod:`orpheus.sn.operators.streaming` continue to wrap their matvec
in :class:`scipy.sparse.linalg.LinearOperator` directly, untouched by
this module. The ORPHEUS↔scipy boundary for the Krylov accelerator is
now internal to :class:`~orpheus.numerics.iteration.KrylovAcceleration`
— a single ravel-aware adapter lifts the flat scipy vector to the typed
carrier, runs the carrier-space matvec, and ravels the result back — so
both the system matvec and the preconditioner route through one site
(single source of truth, #257 S7). Concretely, the within-group system
matvec is the named carrier-space closure ``loss_minus_gains(psi) =
A.apply(psi) - Σ_gains g.apply(psi)`` (it reads like the within-group
operator instead of being buried under ravel plumbing), and *both* it
and the preconditioner are wrapped by the **one** module-private
``_as_scipy_linop(carrier_matvec, template, n)`` adapter, whose
``template`` argument carries the typed carrier so ``from_flat`` /
``to_flat`` round-trip the scipy vector. The flat-only public twin
``as_scipy_linop`` (a zero-production-caller adapter, the degenerate flat
case of the ravel-aware one) was retired in the consolidation; a future
caller needing a standalone operator→scipy adapter (e.g. a DSA
preconditioner built as an operator, #2) should *generalise*
``_as_scipy_linop`` to accept an ``op.apply`` callable, not resurrect the
flat twin.


.. _diagonal-operator:

Diagonal operator on a tagged axis
==================================

The simplest non-trivial operator beyond the composition primitives
is the **diagonal-on-an-axis** operator
:class:`~orpheus.numerics.operator.DiagonalOperator`. For a 1-D
weight vector :math:`w \in \mathbb{R}^N` and target axis ``axis``,
its action on a multi-axis tensor :math:`x` is elementwise
multiplication along ``axis``:

.. (vv-status rationale) Verified by
   ``tests/numerics/test_diagonal_operator.py`` — apply against
   ``np.einsum`` reference on randomised tensors, self-adjointness,
   and round-trip ``apply ↔ solve`` bit-identity.
.. vv-status: diagonal-operator-action documented

.. math::
   :label: diagonal-operator-action

   (D x)_{\ldots,\,n,\,\ldots}
   \;=\; w_n \, x_{\ldots,\,n,\,\ldots}.

All other axes broadcast through unchanged. This is the canonical
"diagonal in some basis" operator — the abstraction Grand Report v3
§9 names :math:`W` (``AngularWeightMatrix``) when the basis is the
discrete-ordinate set of an angular cubature, and the same primitive
any method needs for "multiply-by-weights along one axis":

.. list-table:: Cross-method consumers of DiagonalOperator
   :header-rows: 1
   :widths: 24 38 38

   * - Consumer
     - Role of the diagonal
     - Source
   * - **SN** (:math:`Y^* W` projection)
     - :math:`W` = the angular cubature weights
       :math:`w_n`; the operator multiplies the angular axis of
       :math:`\psi`.
     - This work, Wave 1 (the :math:`W` inside the SH frame's
       analysis face ``frame.analysis``).
   * - **MoC**
     - Track-weight diagonal :math:`w_t` on the track axis of an
       angular flux defined per-track.
     - Future MoC consumer.
   * - **CP**
     - Region-volume diagonal :math:`V_i` on the cell axis of a
       collision-probability matrix.
     - Future CP consumer.
   * - **MC**
     - Importance-weighting diagonal on the source / track axis of
       a tally.
     - Future MC consumer.

Self-adjointness is automatic for real-valued weights:
:meth:`apply_transpose` is the same code path as :meth:`apply`.
Invertibility is by-element: if every weight is non-zero, the
operator reports ``is_invertible = True`` and its :meth:`solve`
divides by :math:`w_n` along ``axis``. If any weight is zero,
``is_invertible`` reports ``False`` and ``inverse()`` / ``solve`` raise
:class:`~orpheus.numerics.operator.NotInvertible` **eagerly** (the
value-dependent arm of Design C) — a zero weight has no inverse, and the
harmful-stub anti-pattern the three-layer surface exists to prevent (a
downstream Krylov consumer silently dividing by zero) is caught upfront.

The implementation does NOT eagerly materialise an
:math:`N \times N` diagonal matrix. The action is a single
broadcast-multiply
(``w.reshape((1, ..., -1, ..., 1)) * x``) so memory cost is
:math:`O(N)` regardless of the input tensor's shape. For the SN
angular axis this matters: an :math:`(N, n_x, n_y, n_g)` field with
:math:`N = O(10^2)` and :math:`n_x \cdot n_y \cdot n_g = O(10^7)`
does NOT need a :math:`(10^7, 10^7)` materialised diagonal.


.. _multiplication-operator-promotion:

The multiplication operator — a coefficient field promoted (§5.7)
=================================================================

The grand report §5.7 closes a loop that the rest of the operator
algebra leaves open. Throughout this page, an operator is a
:class:`~orpheus.numerics.operator.LinearOperator` — an opaque
``Generic[V]`` whose ``apply`` carries a flux to a flux, with nothing in
its type that says *what physics it carries*. For the collision
operator that opacity is a missed opportunity: the collision term
:math:`C\psi = \Sigma_t\,\psi` is **nothing but a cross section acting
by pointwise multiplication**. The cross section is not an input to the
operator; the cross section **is** the operator. §5.7 makes that
identity literal: a :class:`~orpheus.transport.fields.cross_section_field.CrossSectionField`
:math:`f` is *promoted* to the multiplication operator :math:`M[f]`, and
:math:`C = M[\sigma_t]` becomes a named instance of that promotion
rather than an anonymous broadcast-multiply buried in the operator's
``apply`` body (the action now lives once, in
:meth:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator.apply`).

The multiplier-algebra embedding
--------------------------------

The promotion is the **multiplier-algebra embedding** of measure theory:
for the Hilbert space :math:`L^2` of square-integrable flux distributions
over phase space, every bounded measurable coefficient
:math:`f \in L^\infty` defines a bounded operator on :math:`L^2` by
pointwise multiplication,

.. (vv-status rationale) Structural / definitional identity — the
   multiplier-algebra embedding M: L^∞ → B(L²). Not a solver claim; the
   verifiable content is the multiplier-algebra law-suite
   (``tests/transport/test_multiplication_operator.py``) which pins the
   discrete realization :eq:`multiplication-operator-action`.
.. vv-status: multiplication-operator-embedding documented

.. math::
   :label: multiplication-operator-embedding

   M \;:\; L^\infty \;\longrightarrow\; B(L^2),
   \qquad
   (M[f]\,\psi)(\xi) \;=\; f(\xi)\,\psi(\xi),

where :math:`\xi = (\hat\Omega, g, \vec r)` ranges over the discrete
phase space (ordinate, group, cell). The map :math:`M` is the canonical
faithful unital ``*``-homomorphism of :math:`L^\infty` *onto the
diagonal subalgebra* of :math:`B(L^2)` — "diagonal" because
multiplication by :math:`f` couples no two phase-space points: the
output at :math:`\xi` reads the input only *there*. This is the algebraic
content of the §5.6 *locality* criterion that separates the §5.7
multiplication operator from the nonlocal integral kernels
(:ref:`integral-kernel-category`): :math:`M[f]` is the diagonal
sub-algebra, the kernels are everything off-diagonal.

For the leading-ordinate broadcast on the SN per-ordinate carrier
:math:`\psi(\hat\Omega_n, g, \vec r)`, the discrete embedding is the
group-and-space-indexed broadcast over the ordinate axis:

.. (vv-status rationale) #257 S3b — the §5.7 multiplier promotion. The
   broadcast action is verified at the VALUES level against the legacy
   ``σ[None]·ψ`` (0 ULP, ``assert_array_equal`` — the generalized engine
   reduces to the same broadcast-multiply) and the multiplier-algebra
   laws are verified as intrinsic properties (``tests/transport/
   test_multiplication_operator.py``). The structurally-independent
   physics backing is the k_∞ = νΣf/Σa analytical limit and the
   streaming-equilibrium ψ = Q/σ_t reference, both of which route σ_t
   through the promoted C.
.. vv-status: multiplication-operator-action documented

.. math::
   :label: multiplication-operator-action

   (M[f]\,\psi)_{n,g,\vec r} \;=\; f_{g,\vec r}\,\psi_{n,g,\vec r}.

The action is delegated to the N-D
:class:`~orpheus.numerics.operator.DiagonalOperator` broadcast engine
(:eq:`diagonal-operator-action`, #257 S3a) as
``DiagonalOperator(f.values, broadcast_axes=(0,))``: the coefficient is
an :math:`(n_g, *\text{spatial})` field, broadcast over the leading
ordinate axis. The transport-level
:class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`
does **not** re-derive the broadcast — it *builds the engine once* over
its immutable coefficient and forwards (``coding-elegance`` Pattern 2,
single source of truth), so the transport operator and the numerics
engine agree by construction rather than by a copied predicate.

The cross section IS the operator
---------------------------------

The promotion is the moment the opaque ``Generic[V]`` becomes an honest
*carrier*. Before §5.7, the collision operator stored a raw
:math:`\sigma_t` array (an ``ndarray`` with no type-level meaning) and
re-implemented the broadcast in its ``apply`` body — the type
``CollisionOperator`` said nothing, and the meaning lived in a comment.
After §5.7,
:class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`
stores **only** a :class:`~orpheus.transport.fields.cross_section_field.CrossSectionField`
(the cone-typed coefficient of #257 S1), sourced through the typed
:meth:`~orpheus.transport.mesh.material_xs_field.MaterialXSField.total_cross_section_field`
accessor (#257 S2). The operator's *identity* is its coefficient field;
the action follows from the embedding. The collision leaf is now a
**plain** :class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`
carrying :math:`\sigma_t` as its coefficient — :math:`C = M[\sigma_t]` is
literally true, with no SN-specific subtype at all. #261 retired the
former ``CollisionOperator`` thin subclass: once the transport base gained
the optional :attr:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator.space`
(for the W-D composition guard) and the bare-array
:meth:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator.from_mesh`
constructor, the subclass added nothing the base lacked. The ``L + C``
dispatch that assembles the bundled
:class:`~orpheus.sn.operators.streaming.InvertibleOperator` lives **one-directionally**
on :meth:`~orpheus.sn.operators.streaming.StreamingOperator.__add__` (keyed on the
``MultiplicationOperator`` base type): ``L + C`` is the canonical (and
only) ordering, because the ``numerics ↛ transport ↛ sn`` layer order
forbids a transport-level multiplier from dispatching back onto an ``sn``
operator (that would be a ``transport → sn`` upward import).

The multiplier-algebra laws
---------------------------

That :math:`M` is a faithful unital ``*``-homomorphism is not decoration
— it is the set of *intrinsic properties* the promotion must satisfy, and
each is pinned as a law-suite test (the user directive: every
math-bearing type ships a test of its defining laws). The laws, verified
in ``tests/transport/test_multiplication_operator.py`` on a discriminating
``nx=5 ≠ ny=3, ng=2`` carrier with asymmetric heterogeneous coefficients:

.. list-table:: Multiplier-algebra laws (faithful unital ``*``-homomorphism)
   :header-rows: 1
   :widths: 26 30 44

   * - Law
     - Statement
     - Meaning / verification
   * - **Unit**
     - :math:`M[1] = I`
     - The constant-one coefficient is the identity operator (the
       embedding is *unital*).
   * - **Zero (codomain-aware)**
     - :math:`M[0] = 0`
     - The zero coefficient is the zero operator — but its codomain is a
       **source**, not a flux: collision turns flux into a reaction rate,
       so ``M[0].apply`` emits a zeroed
       :class:`~orpheus.transport.source_sinks.angular_source_sink.AngularSourceSink`.
   * - **Linearity**
     - :math:`M[af + bg] = a\,M[f] + b\,M[g]`
     - :math:`M` is linear in the coefficient (the embedding is a vector-
       space map). Tested on ≥2-group asymmetric heterogeneous fields.
   * - **Homomorphism**
     - :math:`M[f]\,M[g] = M[f g]`
     - Composing two multiplications multiplies the coefficients. Tested
       at the VALUES level on the raw :math:`\sigma\cdot\sigma'` product
       (which has units :math:`\mathrm{cm}^{-2}` — the units grading is
       exactly why coefficient-field ``*`` is deferred to the values
       layer rather than the typed field layer).
   * - **Self-adjointness**
     - :math:`M[f]^* = M[f]` for real :math:`f`
     - A real-valued multiplication is self-adjoint;
       ``apply_transpose`` is the same code path as ``apply``.
   * - **Spectrum**
     - :math:`\mathrm{spec}(M[f]) = \mathrm{ess\,range}(f)`
     - :math:`M[f]` is invertible **iff** :math:`f` is bounded away from
       zero; the inverse is :math:`M[1/f]`. This is the honest
       invertibility gate (below).

The spectrum law and the honest invertibility gate
--------------------------------------------------

The spectrum of a multiplication operator is the essential range of its
coefficient: :math:`M[f]` has a bounded inverse iff
:math:`\mathrm{ess\,inf}\,|f| > 0`, in which case
:math:`M[f]^{-1} = M[1/f]`. The promotion enforces this at the value
level —
:attr:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator.is_invertible`
reports ``True`` **iff** :math:`\min|f| > 0`, single-sourced from the
broadcast engine (the operator reads the engine's ``is_invertible``, so
the transport operator and the numerics engine agree by construction).
The multiplier is the *value-dependent* arm of Design C
(:ref:`design-c-structural-value-split`): it always **declares**
``inverse()`` / ``solve`` (the type *can* invert), and refuses eagerly.

This is a **behavioral strengthening**, an illegal-states-unrepresentable
hardening (``coding-elegance`` Pattern 4). The legacy ``CollisionOperator``
advertised a ``solve``
*unconditionally*; on a :math:`\sigma = 0` entry its ``solve`` divided
:math:`q/\sigma` and produced a **silent IEEE NaN** that propagated into
the iterate. The promotion refuses the inverse it does not have: a zero
coefficient has no inverse, so
:attr:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator.is_invertible`
reports ``False`` and any request for the inverse
(``inverse()`` / ``solve``) raises
:class:`~orpheus.numerics.operator.NotInvertible` **eagerly** rather than
producing a NaN at *call* time. Construction still succeeds (the gate
governs only the inverse, never blocks the object); ``apply`` is
unaffected.

.. note::

   The gate change is purely additive honesty — an audit (#257 S3b)
   confirmed **no production path** relies on a :math:`\sigma=0` collision
   ``solve``. Since B.2d d3 the within-group :math:`L + C` composition is
   spelled in **one** place — the fused-factor primitive
   :func:`~orpheus.sn.coupled_system.build_streaming_collision`, called by
   BOTH :func:`~orpheus.sn.coupled_system.build_within_group_system` (whose
   :class:`~orpheus.sn.coupled_system.WithinGroupSystem` record the
   solver's within-group inners consume) and the solver's own ``L``
   binding legs (the former per-solver ``Streaming + Collision``
   spellings — the third copy was the L-002 collapse trigger). That builder
   sources
   :math:`\sigma_t > 0` (total cross sections are bounded away from zero),
   and the removal cross section :math:`\sigma_r` ``solve`` appears only in
   a docstring, with no live caller. The bundled
   :class:`~orpheus.sn.operators.streaming.InvertibleOperator` has its own
   stricter construction-time ``min σ > 0`` guard, consistent with the new
   gate.


.. _functional-category:

The functional category (the §5.6 suffix law)
=============================================

The grand report §5.6 *suffix law* partitions the maps of the algebra
by what they return: an **Operator** maps a field to a field
(:eq:`operator-apply`), a **Kernel** integrates a field against a
measure (:eq:`scattering-as-tensor-product-sum`), and a **Functional**
maps a field to a **scalar** — or, fiberwise over space, to a
scalar-*field*. The functional surface is the single method
``evaluate(x) -> R``; it deliberately carries **none** of the
:class:`~orpheus.numerics.operator.LinearOperator` surface (no
``apply``, no ``is_invertible`` / ``is_adjointable``), and that
disjointness is the category's defining property (#257 S5,
:class:`orpheus.numerics.functional.Functional`).

The category is seated at **both layers**. The generic numerics floor is
:class:`~orpheus.numerics.functional.InnerProductFunctional`\ ``(weight,
axis)`` — the co-vector :math:`\langle w, \cdot\rangle` whose ``evaluate``
is the single primitive ``(w * x).sum(axis, keepdims=True)``. The
transport leaf
:class:`~orpheus.transport.reaction_rate_functional.ReactionRateFunctional`
**specialises** it (carrying the weight as a domain-typed
:class:`~orpheus.transport.fields.cross_section_field.CrossSectionField`,
which brings ``.mesh`` and the ``1/cm`` units), exactly as
:class:`~orpheus.transport.frames.harmonic_frame.HarmonicFrame` specialises
:class:`~orpheus.numerics.frame.GalerkinFrame`. The canonical instance is
the per-cell reaction-rate **density**, contracting a reaction cross
section against the flux over the source groups:

.. math::
   :label: production-rate-functional

   r_x(\vec r) \;=\; \langle \Sigma_x, \phi\rangle
   \;=\; \sum_{g'} \Sigma_{x,g'}(\vec r)\,\phi_{g'}(\vec r),

with the group axis collapsed and **no** volume measure folded in (it
is the density, not the integral :math:`\int r_x\,dV`). The two named
instances are the **production rate** (:math:`\Sigma_x = \nu\Sigma_f`)
and the **absorption rate** (:math:`\Sigma_x = \Sigma_a`); the
Rayleigh-quotient eigenvalue is their ratio
:math:`k = \langle\nu\Sigma_f,\phi\rangle / \langle\Sigma_a,\phi\rangle`.

This functional is not a parallel *description* of a contraction the
operator algebra performs elsewhere — it **is** the contraction. It is
the row-factor :math:`\langle\nu\Sigma_f|` of the fission rank-1 dyad
:math:`F = |\chi\rangle\langle\nu\Sigma_f|`
(:ref:`fission-as-dyad`): :meth:`RankOneOperator.apply
<orpheus.numerics.operator.RankOneOperator.apply>` routes the fission
matvec *through* ``functional.evaluate``, so there is no separate "fused"
realization to drift from the named factor (the procedural twin the
earlier S5 design carried is **dissolved** — :ref:`fission-as-dyad`
records the upgrade). Naming the contraction turns the most physically
central diagnostic in criticality into a typed, inspectable Functional.

The volume-integrated companion
:class:`~orpheus.transport.reaction_rate_functional.IntegratedReactionRate`
returns the **scalar** :math:`R_x = \int_V \langle\Sigma_x,\phi\rangle\,dV`
— the per-cell density above, integrated against the mesh's canonical
``volume_measure`` (single source: it reuses
:eq:`production-rate-functional`'s density, no independent re-derivation
of either reduction). It is the typed object behind the
:math:`k`-eigenvalue numerator and denominator
(:ref:`integrated-reaction-rate-keff`), and the :math:`\phi^\dagger=1`
degenerate of the homogenisation Petrov–Galerkin campaign's
adjoint-weighted bilinear :math:`\langle\phi^\dagger, M[\Sigma_x]\,\phi\rangle`
(a future adjoint flux replaces the implicit :math:`\phi^\dagger = 1`).

.. vv-status: functional-category documented
.. vv-status: production-rate-functional documented

The §5.6 suffix law as a type-system fact
-----------------------------------------

The §5.6 *suffix law* is not a taxonomy imposed on the prose — it is a
structural fact the type system enforces by what each category's surface
*is*. Three categories partition the maps of the algebra by their
codomain:

.. list-table:: The §5.6 suffix law (the codomain partition)
   :header-rows: 1
   :widths: 18 24 24 34

   * - Category
     - Signature
     - Type relationship
     - ORPHEUS realization
   * - **Operator**
     - field :math:`\to` field
     - the base —
       :class:`~orpheus.numerics.operator.LinearOperator`
     - streaming :math:`L`, collision :math:`C = M[\sigma_t]`
       (:eq:`multiplication-operator-action`)
   * - **Kernel**
     - field :math:`\to` field, *nonlocal*
     - a **refinement** of Operator (adds ``kernel``)
     - scattering :math:`S`, fission :math:`F`
       (:ref:`integral-kernel-category`)
   * - **Functional**
     - field :math:`\to` scalar
     - a **disjoint sibling** of Operator (shares no member)
     - the reaction-rate density :math:`r_x(\vec r)`
       (:eq:`production-rate-functional`)

The asymmetry between the two non-Operator categories is the point. A
**Kernel** *is* a :class:`~orpheus.numerics.operator.LinearOperator` (it
still maps a field to a field, it still has ``apply`` and the two-axis
predicates) and merely *adds* the ``kernel`` member — so it is a
strict refinement (:ref:`integral-kernel-category`). A **Functional** is
*not* a :class:`~orpheus.numerics.operator.LinearOperator` at all: its
sole surface is ``evaluate(x) -> R``, and it deliberately carries
**none** of the operator surface — no ``apply``, no ``is_invertible`` /
``is_adjointable``, no ``solve``, no ``apply_transpose``, no ``.H``, no
``domain`` / ``codomain``.
That *disjointness* is the category's defining property, not an
omission, and it is checkable: ``isinstance(p, Functional)`` is true while
``isinstance(p, LinearOperator)`` is false, and the discriminator foils
both directions (a bare operator is not a Functional; the production-rate
functional is not an operator). The category is the
:class:`orpheus.numerics.functional.Functional` Protocol (#257 S5, the L1
numerics floor — the **co-vector companion of**
:class:`~orpheus.numerics.vector.Vector`: a ``Vector`` is what an
operator acts *on*, a ``Functional`` is what acts *on* a vector to
produce a scalar).

Variance: a Functional is contravariant in, covariant out
---------------------------------------------------------

The Functional's typevars carry a *different* variance profile than the
operator's, and the difference is structural. An operator's
``apply(x: V) -> V`` uses :math:`V` **both ways** (it consumes a flux and
produces a flux), so :math:`V` is *invariant by dual use* — and pyright
emits no variance warning. A functional's ``evaluate(x: V_contra) ->
R_co`` uses its input **only** as a consumer (contravariant) and its
result **only** as a producer (covariant). Declaring the Functional over
the operator's invariant :math:`V` would be a type error of variance;
the correct declaration is a contravariant input typevar and a
covariant, *unbounded* result typevar (the result is unbounded because
the production rate's "scalar" is fiberwise a scalar-**field** over space,
not a Python :class:`float` — a ``float | V`` union would mistype it).
This is the co-vector mirror of the invariant :math:`V` and is the
load-bearing typing lesson of S5: the variance is not cosmetic, it is the
formal statement that a functional is a co-vector.

The reaction rate is a density, not an integral
-----------------------------------------------

The :class:`~orpheus.transport.reaction_rate_functional.ReactionRateFunctional`
of :eq:`production-rate-functional` returns the per-cell reaction-rate
**density** :math:`r_x(\vec r)` — the cross section contracted against the
flux over the source groups, group axis collapsed, **no volume measure
folded in**. The "no measure" choice is a deliberate Mode-3
(missing-factor) guard turned into a *named contract*: :math:`r_x` is the
density :math:`\Sigma_x\,\phi`, not the integral :math:`\int r_x\,dV`. A
consumer that needs the *integrated* rate uses the dedicated
:class:`~orpheus.transport.reaction_rate_functional.IntegratedReactionRate`
(which applies the mesh ``volume_measure`` to this exact density); folding
the measure into :math:`r_x` itself would silently double-count it the
moment a second consumer integrated again. The split is the structural
division between the two functionals — the density carries the group
contraction, the integrated rate adds the spatial measure on top of it.
The coefficient side is the cone-typed
:class:`~orpheus.transport.fields.cross_section_field.CrossSectionField`
(#257 S1) carrying :math:`\Sigma_x` through the typed
:meth:`~orpheus.transport.mesh.material_xs_field.MaterialXSField.fission_production_field`
accessor (#257 S2) for production, or the absorption accessor for
:math:`\Sigma_a`.

.. _integrated-reaction-rate-keff:

The eigenvalue numerator and denominator are integrated reaction rates
----------------------------------------------------------------------

The SN :math:`k`-eigenvalue is routed through
:class:`~orpheus.transport.reaction_rate_functional.IntegratedReactionRate`
directly:

.. math::
   :label: keff-as-integrated-rates

   k \;=\; \frac{\displaystyle\int_V \langle\nu\Sigma_f,\phi\rangle\,dV
                  \;+\; (n,2n)}
                {\displaystyle\int_V \langle\Sigma_a,\phi\rangle\,dV} ,

with both the numerator's fission term and the denominator the
volume-integrated reaction rate :math:`R_x = \int_V\langle\Sigma_x,\phi
\rangle\,dV` of :eq:`production-rate-functional`
(:meth:`SNSolver.compute_keff <orpheus.sn.solver.SNSolver.compute_keff>`
over :meth:`~orpheus.sn.solver.SNSolver.compute_production_rate`). The
:math:`(n,2n)` channel is an **explicit additive term** — a second
neutron-multiplying reaction, *not* a :math:`\langle\Sigma_x,\phi\rangle`
rate, so it is added on top rather than folded into a cross section — and
is exactly zero on a no-:math:`(n,2n)` mixture. Reusing
:class:`IntegratedReactionRate
<orpheus.transport.reaction_rate_functional.IntegratedReactionRate>` for
both production and absorption gives the eigenvalue numerator and
denominator a single typed source for the :math:`\langle\Sigma_x,\phi
\rangle` contraction and its volume integral, the
:math:`\phi^\dagger\!=\!1` degenerate of the homogenisation
Petrov–Galerkin bilinear :math:`\langle\phi^\dagger, M[\Sigma_x]\,\phi
\rangle`.

.. note::

   **Scope.** The SN :math:`k`-eigenvalue routes through
   :class:`IntegratedReactionRate
   <orpheus.transport.reaction_rate_functional.IntegratedReactionRate>`
   for **both** its numerator and denominator. Diffusion (#290) routes
   its **production** rate through the same functional (the #270
   diffusion arm — see :doc:`/theory/diffusion_1d`), but poses its
   denominator as the integrated loss-operator action
   :math:`\langle 1, (A\psi)_{\rm bulk}\rangle_V` (= absorption +
   leakage by the column-sum theorem), not a second reaction-rate
   contraction. The CP and MoC eigenvalues are **not yet** routed;
   the homogeneous 0-D case has no spatial integral to fold, so it
   does not participate. Do not read this section as a claim that every
   solver's :math:`k` flows through the typed functional in the same
   way — it is the SN path, with diffusion's production numerator the
   one other current consumer.

.. (vv-status rationale) Structural / decomposition label: the SN
   eigenvalue expressed as a ratio of integrated reaction rates. Not a new
   solver claim (the SN keff value is verified by the existing eigenvalue
   gates); the verifiable content of the *routing* is the closed-form
   k∞-as-ratio gate ``tests/transport/test_integrated_reaction_rate.py``.
.. vv-status: keff-as-integrated-rates documented

.. _reaction-rate-kinf-oracle-section:

Verification: the closed-form :math:`k_\infty` oracle, not a procedural twin
----------------------------------------------------------------------------

The earlier S5 design verified the functional **byte-for-byte (0 ULP)
against the rank-1 ``inner`` reduction** :math:`\chi\cdot\mathrm{evaluate}
\equiv F.\mathrm{apply}`. That cross-check was *procedurally* independent
(two code paths) but **not** *structurally* independent (both sides ran
the same NumPy ``(w * x).sum(axis, keepdims)`` primitive on the same
arrays) — by ``vv-principles`` L11 a procedural twin proves the two
spellings agree, not that either is **correct**. With the matvec now
routed *through* ``functional.evaluate`` (:ref:`fission-as-dyad`) the twin
no longer even exists to compare against: there is one contraction, not
two.

The correctness floor is therefore re-anchored on a **structurally
independent closed-form ground** — the infinite-medium decomposition

.. math::
   :label: reaction-rate-kinf-oracle

   k_\infty \;=\; \lambda_{\max}\!\bigl(A^{-1}F\bigr),
   \qquad A = \mathrm{diag}(\Sigma_a) - \Sigma_{s}, \quad
   F = |\chi\rangle\langle\nu\Sigma_f| ,

whose dominant eigenpair :math:`(k_\infty, \phi^*)` comes from a
:func:`numpy.linalg.eig` of the transfer matrix — a path that shares **no
primitive** with the rank-1 functional. The
:class:`~orpheus.transport.reaction_rate_functional.ReactionRateFunctional`
is pinned by evaluating production **and** absorption *independently* at
the converged spectrum :math:`\phi^*` and checking
:math:`\langle\nu\Sigma_f,\phi^*\rangle / \langle\Sigma_a,\phi^*\rangle =
k_\infty` (``tests/transport/test_reaction_rate_functional.py``). Pinning
the two functionals *separately* (not just their ratio) is what gives the
gate teeth a ratio test lacks: a shared-factor error that scales both the
numerator and the denominator — a mis-scaled accessor, a spurious volume
fold on both — cancels in :math:`k` but is caught term-by-term.

This is a legal **eigenvalue claim** paired with the **closed-form**
pillar (``vv-principles``: eigenvalue claims need a closed-form or
semi-analytical reference, never MMS). It carries genuine **flux-shape
teeth** because the test uses a **4-group** mixture whose converged
:math:`\phi^*` is genuinely non-flat — the 2-group case is degenerate (its
:math:`\phi^*` is coincidentally flat ``[0.707, 0.707]``, so its flat-flux
ratio equals :math:`k_\infty` and is flux-shape-blind), so the 4-group leg
is **mandatory** (``vv-principles`` anti-pattern #3: a 1-group, or any
flat-spectrum, eigenvalue test cannot detect a flux-shape error). A
second leg pins the per-cell density against an explicit Python
double-loop (``hand_derived_production_density`` — an accumulation with no
NumPy reduction, the L11 structurally-independent reference), and a third
asserts the **no-volume-measure** contract (Mode-3). Together they replace
the retired procedural twin with a reference at the right level for each
claim.

.. (vv-status rationale) Structural / decomposition label: the closed-form
   k∞ = λ_max(A⁻¹F) identity that grounds the reaction-rate functional's
   correctness. Not itself a solver claim; the verifiable content is the
   eigenvalue/closed-form gate ``tests/transport/test_reaction_rate_functional.py``
   (production AND absorption pinned independently at φ*) and the
   ``tests/transport/test_integrated_reaction_rate.py`` k∞-as-ratio gate.
.. vv-status: reaction-rate-kinf-oracle documented

Why the estimators are NOT Functionals
---------------------------------------

The criticality eigenvalue and production-rate **estimators**
(:meth:`~orpheus.numerics.iteration.KEigenvalue.compute_keff` and
:meth:`~orpheus.numerics.iteration.KEigenvalue.compute_production_rate`)
are the obvious candidates to wrap as ``Functional`` objects — and they
are deliberately **not**. The eigenvalue estimator is a **ratio** of two
triple-dependent contractions,
:math:`\sum(F\psi)\,/\,(\sum(A\psi) - \sum(S\psi))` — it consumes the
whole operator triple :math:`(A, S, F)` (carried on the ``KEigenvalue``
instance) together with the iterate :math:`\psi`, not a lone field acted
on by a single co-vector. That ratio-of-triple-contractions shape is not
the ``evaluate(x) -> R`` shape of a
:class:`~orpheus.numerics.functional.Functional` (a linear co-vector on
one vector space). The category simply now *names* what their
field-to-scalar **core** is (the production-rate contraction
:math:`\sum(F\psi)`), without forcing the estimators into a wrapper that
would misrepresent their arity. They stay bare (hardwired) methods,
arithmetic bit-identical to the pre-R8 module-level defaults they
replaced (pinned by ``tests/numerics/test_estimators_as_functionals.py``),
and the honesty of *not* wrapping them is itself a category-correctness
claim.

.. note:: **Injection seam retired (R8, #259 P1, 2026-07-03).** Before
   this, these estimators were injectable *callables* — the
   ``_default_keff_estimator`` / ``_default_production_estimator`` module
   functions were the **defaults** of ``KEigenvalue``'s ``keff_estimator``
   / ``production_estimator`` kwargs, which a caller could override.  The
   kwargs, the ``KeffEstimator`` / ``ProductionEstimator`` aliases, and
   the ``_default_*`` functions are gone; the spellings moved verbatim
   onto the hardwired methods above.  The seam was dead **by design** —
   at a converged eigenpair every estimator consistent with the posed
   problem agrees, so an injected *different* estimator could only be
   inconsistent (illegal states unrepresentable).  See the full R8
   rationale — the consistency theorem and the honest-``A.apply``
   contract — in the *KEigenvalue: outer power iteration* section of
   :doc:`discrete_ordinates`.

.. note::

   **Source map.** Category Protocol:
   :class:`orpheus.numerics.functional.Functional` (L1). Generic numerics
   leaf: :class:`~orpheus.numerics.functional.InnerProductFunctional`
   ``⟨w, ·⟩`` (L1). Transport leaves:
   :class:`~orpheus.transport.reaction_rate_functional.ReactionRateFunctional`
   (the per-cell density, specialising ``InnerProductFunctional``) and
   :class:`~orpheus.transport.reaction_rate_functional.IntegratedReactionRate`
   (the volume-integrated scalar) (L2). Rank-1 constructor:
   :func:`~orpheus.numerics.operator.outer`. Intrinsic-category gate:
   ``tests/transport/test_functional_category.py`` (Functional ≠
   LinearOperator, both directions + discriminator foils). Correctness:
   ``tests/transport/test_reaction_rate_functional.py`` (the closed-form
   :math:`k_\infty = \lambda_{\max}(A^{-1}F)` per-term oracle —
   production AND absorption pinned independently — + the hand-derived
   double-loop density reference + the no-measure guard);
   ``tests/transport/test_integrated_reaction_rate.py``
   (:math:`k_\infty` as the ratio of integrated rates, incl. the
   ``(n,2n)`` term). Dyad laws: ``tests/numerics/test_outer_dyad.py``
   (action / rank-1 / adjointability / linearity). Estimator honesty:
   ``tests/numerics/test_estimators_as_functionals.py``. The
   field-to-scalar contraction is coded in three places today
   (SN / CP / numerics) — unifying that fragmentation is tracked on #259.


.. _integral-kernel-category:

The integral-kernel category (the §5.6 Kernel suffix)
=====================================================

The §5.6 suffix law's middle term is the **Kernel**: a
:class:`~orpheus.numerics.operator.LinearOperator` whose action is
**nonlocal** — it integrates the carrier field against a measure on one
or more axes — distinct from a LOCAL / diagonal
:class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`
(:eq:`multiplication-operator-action`, the §5.7 Operator) and from a
field→scalar :class:`~orpheus.numerics.functional.Functional`
(:eq:`production-rate-functional`). The single discriminator is
**locality** (Frame 3): a multiplication operator's output at a point is
a pointwise function of the input *there*, while a Kernel's output reads
the input across an integrated axis. #257 S6 names this category as the
:class:`orpheus.transport.operators.integral_kernel_operator.IntegralKernelOperator`
Protocol — a **refinement of** LinearOperator (it still has ``apply`` +
the two-axis predicates, UNLIKE the disjoint Functional) that adds a
single ``kernel`` member exposing the integral structure as a
:class:`~orpheus.numerics.operator.LinearOperator`:

.. math::
   :label: integral-kernel-category

   (A\,\psi)(x) \;=\; \int K(x, x')\,\psi(x')\,d\mu(x') ,
   \qquad K \;=\; A.\mathrm{kernel}.

The two named transport instances are the Boltzmann emission kernels.
**Fission** exposes the rank-1 **dyad** :math:`F = |\chi\rangle\langle
\nu\Sigma_f| = \texttt{outer}(\chi,\,\mathrm{ReactionRateFunctional}(\nu
\Sigma_f))` (a :class:`~orpheus.numerics.operator.RankOneOperator`, lifted
to a :class:`~orpheus.numerics.operator.TensorProductOperator` only to
advertise the spatial-axis broadcast). Its row co-vector
:math:`\langle\nu\Sigma_f|` is the S5 :eq:`production-rate-functional`
density, and the matvec routes *through* that functional's ``evaluate``
(:ref:`fission-as-dyad`) — there is no fused procedural twin. **Scattering**
exposes the anisotropic Legendre redistribution :math:`R \circ \Lambda
\circ M` (an :class:`~orpheus.numerics.operator.OperatorProduct`, the
genuinely-nonlocal-in-angle part of
:eq:`scattering-as-tensor-product-sum`); the isotropic :math:`P_0`
in-scatter and the :math:`(n,2n)` doubling are the local / separate
components of the full scattering ``apply``. Fission is the **rank-1
(single-mode) degenerate** of scattering's multi-mode spectral sum — the
polyadic/block-term view of :ref:`emission-kernels-btd` makes the
relationship precise. The matvec arms of both operators are UNCHANGED in
S6 (additive, bit-identical).

The scattering kernel :math:`R\circ\Lambda\circ M` is not an arbitrary
nonlocal operator: it is the **spectral theorem** :math:`A = U\Sigma U^*`
of a rotationally-invariant kernel. The scattering kernel
:math:`\Sigma_s(\hat\Omega\cdot\hat\Omega')` depends only on the
direction cosine — a *zonal* kernel — so by the Funk–Hecke theorem the
spherical harmonics are its eigenfunctions, with eigenvalues the
Legendre moments :math:`\Sigma_{s,\ell}` (the diagonal of
:math:`\Lambda` =
:class:`~orpheus.transport.operators.scattering.LegendreMomentScattering`). Reading
:math:`M` = :math:`U^*` (change of basis *into* the eigenbasis),
:math:`\Lambda` = :math:`\Sigma` (the diagonal spectrum), and
:math:`R` = :math:`U` (synthesis *out of* it) is what makes the
conjugation
:math:`S = \tfrac{1}{W}\,\texttt{frame.conjugate}(\Lambda)` — the
scattering **2-cell** of the carrier-grid double category
(:ref:`carrier-grid-double-category`) — a *spectral* statement: the
horizontal adjoint pair :math:`(M, R)` is the unitary diagonalising the
vertical :math:`\Lambda`. This is also why the frame is **owned by**
the scattering operator (its constructor binds the frame order to
``scattering_order``): the spherical harmonics are *scattering's*
eigenbasis, and the streaming operator — the :math:`\ell=1` direction
irrep — is the one transport operator the basis does **not**
diagonalise (it couples :math:`\ell\!\leftrightarrow\!\ell\pm1`, the Pℓ
recurrence). The full Funk–Hecke / Schur derivation, the literature
corroboration, and the unifying principle *"an operator owns its frame
iff the frame is its eigenbasis"* (which also explains why energy
condensation and spatial homogenisation are Petrov-Galerkin, not
Galerkin) are in :ref:`frame-eigenbasis-ownership`
(:doc:`galerkin_projection`).

.. vv-status: integral-kernel-category documented

The locality criterion completes the partition
-----------------------------------------------

With the Kernel category named, the §5.6 suffix-law partition is
**complete and exhaustive**: every map in the transport algebra is an
Operator, a Kernel, or a Functional, and the boundary between Operator
and Kernel is a single mathematical criterion — **locality**. A
multiplication operator (:eq:`multiplication-operator-action`) is the
**diagonal** sub-algebra: its output at a phase-space point reads the
input only *there*. An integral kernel is everything **off-diagonal**:
its output at :math:`x` integrates the input across an axis,
:math:`(A\psi)(x) = \int K(x,x')\,\psi(x')\,d\mu(x')`. The Kernel's
``kernel`` member is exactly that integrating object :math:`K`, exposed
as a :class:`~orpheus.numerics.operator.LinearOperator` in its own right.

This is why the category is a **refinement** of
:class:`~orpheus.numerics.operator.LinearOperator` and not a disjoint
sibling like the Functional. The Protocol is written
``IntegralKernelOperator(LinearOperator[V], Protocol[V])`` — a Kernel
*is-a* LinearOperator (it keeps the inherited surface: ``apply``, the
``is_invertible`` / ``is_adjointable`` predicates, ``domain`` /
``codomain``) and merely *adds* the ``kernel`` property, which the
``@runtime_checkable`` ``isinstance`` sees on top of the inherited
members. The
refinement is **strict**, which the intrinsic gate verifies in both
directions: a kernel-less LinearOperator
(:class:`~orpheus.numerics.operator.IdentityOperator`, the local
:class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`)
is **not** a Kernel, and a :class:`~orpheus.numerics.functional.Functional`
(no ``apply``, no ``kernel``) is **not** a Kernel either. Variance is
inherited verbatim from LinearOperator — a Kernel's ``apply(x: V) -> V``
is invariant by dual use (no variance warning), *unlike* the co-vector
Functional that needed contravariant input + covariant output.

The two emission kernels: fission and scattering
------------------------------------------------

The two named transport instances are the Boltzmann **emission**
kernels, living in :mod:`orpheus.transport.operators`. Each gains a
``kernel`` member that satisfies the Protocol. **Scattering** is reframed
*in place* and additively: its ``kernel`` is the semantic reading of an
existing matvec that stays byte-for-byte unchanged, a cross-check rather
than a rewrite. **Fission**, by contrast, was genuinely *re-realized* as
the dyad :math:`|\chi\rangle\langle\nu\Sigma_f|` whose ``apply`` routes
through the production-rate functional (:ref:`fission-as-dyad`) — the
earlier "fused realization vs. semantic decomposition" split is gone,
because the dyad *is* the realization. The matvec **value** is unchanged
(0 ULP, :ref:`reaction-rate-kinf-oracle-section`); the realization is no
longer a procedural twin of a separately-named decomposition.

.. list-table:: The two §5.6 emission kernels
   :header-rows: 1
   :widths: 16 38 46

   * - Operator
     - ``kernel`` shape
     - Why this shape
   * - **Fission** :math:`F`
     - rank-1 dyad
       :math:`|\chi\rangle\langle\nu\Sigma_f|`
       (:class:`~orpheus.numerics.operator.RankOneOperator` via
       :func:`~orpheus.numerics.operator.outer`, lifted to a
       :class:`~orpheus.numerics.operator.TensorProductOperator` for the
       spatial-axis broadcast)
     - Fission emits an isotropic source whose group spectrum
       :math:`\chi` is the same everywhere fission occurs — a rank-1
       outer product of the emission-spectrum column with the
       production-rate row co-vector.
   * - **Scattering** :math:`S`
     - :math:`R \circ \Lambda \circ M`
       (an :class:`~orpheus.numerics.operator.OperatorProduct`)
     - The anisotropic :math:`\ell \ge 1` Legendre redistribution: project
       the angular flux onto moments (:math:`M`), apply the per-:math:`\ell`
       group-to-group transfer (:math:`\Lambda`), reconstruct to
       ordinates (:math:`R`). Genuinely nonlocal in angle.

.. _fission-as-dyad:

Fission as the rank-1 dyad :math:`|\chi\rangle\langle\nu\Sigma_f|`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fission is the rank-1 dyad — a reconstruction **column** :math:`|\chi\rangle`
(the emission spectrum) tensored with a functional **row**
:math:`\langle\nu\Sigma_f|` (the production-rate co-vector):

.. (vv-status rationale) Structural / decomposition of the fission
   operator — the rank-1 dyad reading. Not a solver claim; the verifiable
   content is the closed-form k∞ = λ_max(A⁻¹F) correctness of the row
   co-vector (``tests/transport/test_reaction_rate_functional.py``) and the
   0-ULP equivalence of the dyad apply to the matvec arm
   (``tests/sn/operators/test_fission_kernel_crosscheck.py``).
.. vv-status: fission-as-dyad documented

.. math::
   :label: fission-as-dyad

   F \;=\; |\chi\rangle\langle\nu\Sigma_f|
     \;=\; \texttt{outer}\bigl(\chi,\;
           \mathrm{ReactionRateFunctional}(\nu\Sigma_f)\bigr),

read right-to-left as the dyad action :math:`F\,\phi = \chi \cdot
\langle\nu\Sigma_f,\phi\rangle`: the row co-vector
:math:`\langle\nu\Sigma_f|` contracts the flux over groups to the scalar
emission **density** (the S5 :eq:`production-rate-functional`, exposed as
:attr:`FissionOperator.production_rate <orpheus.transport.operators.fission.FissionOperator.production_rate>`),
and the column :math:`\chi` broadcasts it back across the emission groups.
The realization is literally that dyad — ``outer(self.chi,
self.production_rate) & IdentityOperator()`` — and the matvec **routes
through** the row's ``evaluate`` (``coding-elegance`` Pattern 5: build the
*right primitive*; here the right primitive *is* the contraction, not a
parallel description of one). Because there is one contraction and not
two, the "fused vs. unfolded" distinction the earlier S5 design carried
**no longer exists** — :ref:`the verification section <reaction-rate-kinf-oracle-section>`
records the resulting upgrade from the dissolved procedural twin to the
closed-form :math:`k_\infty = \lambda_{\max}(A^{-1}F)` oracle.

.. note::

   **What was tried and rejected — the three-factor composition** :math:`F
   = M_\chi \circ \mathrm{ProductionRate} \circ M_{\nu\Sigma_f}`. An
   earlier reading framed fission as a literal three-operator product:
   multiply by :math:`\nu\Sigma_f`, contract to the density, re-broadcast
   through :math:`\chi`. It was abandoned for two reasons. First, it is
   **dimensionally inconsistent as a composition of operators**: the
   middle factor :math:`\mathrm{ProductionRate}` is a *Functional* (field
   :math:`\to` scalar-field), not a field-:math:`\to`-field operator, so
   the three pieces do not chain as a clean
   :class:`~orpheus.numerics.operator.OperatorProduct` of
   :class:`~orpheus.numerics.operator.LinearOperator`\ s — the
   "composition" silently changed category mid-chain. Second,
   :math:`M_\chi` is **rank-changing** (it takes one scalar density per
   cell and expands it across the :math:`n_g` emission groups), not a
   diagonal
   :class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`,
   so even reading it generously it is not the multiplication operator the
   :math:`M_\bullet` notation implies. The honest form is the **two-factor
   dyad** :math:`|\chi\rangle\langle\nu\Sigma_f|`: the column and the row
   are the only two objects, and the rank-1
   :class:`~orpheus.numerics.operator.RankOneOperator` is the principled
   primitive that *is* the dyad's effect — its ``apply`` routes the matvec
   through the row functional, with no intermediate operator stages to
   materialise.

Scattering as the nonlocal-in-angle kernel
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The scattering :attr:`kernel <orpheus.transport.operators.scattering.ScatteringOperator.kernel>`
is the :class:`~orpheus.numerics.operator.OperatorProduct`
:math:`R \circ \Lambda \circ M` built from the SH frame's analysis
face ``frame.analysis`` (:math:`M`),
:class:`~orpheus.transport.operators.scattering.LegendreMomentScattering` (:math:`\Lambda`),
and the frame's reconstruction face ``frame.reconstruction``
(:math:`R`). It is the strictly-anisotropic :math:`\ell \ge 1` part of the
full scattering ``apply``: the isotropic :math:`P_0` in-scatter, the
:math:`(n,2n)` doubling, and the per-ordinate :math:`1/W` normalisation
are the **local / separate** components that live *outside* the kernel
(a strict sub-component, pinned). The kernel reproduces the existing
anisotropic moment path :math:`R(\Lambda(M\psi))` byte-for-byte (0 ULP);
its physics L1 backing is the existing anisotropic MMS gate
``tests/sn/verification/mms/test_curvilinear_aniso_scattering_p1.py``,
not a new reference.

.. note::

   The full Pℓ form :math:`S = \sum_\ell P_\ell \otimes \Sigma_{s,\ell}`
   (:eq:`scattering-as-tensor-product-sum`) was **considered and
   rejected** as the kernel shape. A
   :class:`~orpheus.numerics.operator.SumOfTensorProductsOperator` would
   require :math:`R` and :math:`M` to be tensor-product *factors* on
   independent axes, but they are rank-changing einsums that mix the
   ordinate and harmonic-coefficient axes — they are *not* valid
   tensor-product factors. Expressing the kernel as a
   ``SumOfTensorProductsOperator`` would be a re-derivation of the
   moment redistribution, not a re-presentation of the existing one;
   that path is tracked as the un-orphaning of
   :class:`~orpheus.numerics.operator.SumOfTensorProductsOperator` on
   #260. The :math:`1/W` per-ordinate normalisation lives *outside* the
   kernel (the kernel is the redistribution, not the quadrature
   weighting).

.. _emission-kernels-btd:

The two emission kernels as a tensor decomposition (a lens, not a type)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Seen through the vocabulary of tensor decompositions, fission and
scattering are the **same object at two ranks**, and naming the relationship
sharpens why fission is the simpler one. This subsection is a *reading* —
a lens that organises the two kernels — **not** a new type to build. The
right-grained primitives already exist (the dyad
:func:`~orpheus.numerics.operator.outer`, the
:class:`~orpheus.transport.frames.harmonic_frame.HarmonicFrame`, and the
``⊗`` / :class:`~orpheus.numerics.operator.TensorProductOperator`); the
decomposition vocabulary names what they already *are*.

.. list-table:: The emission kernels in tensor-decomposition vocabulary
   :header-rows: 1
   :widths: 18 30 52

   * - Kernel
     - Decomposition
     - Reading
   * - **Fission** :math:`F`
     - rank-1 (CP atom)
       :math:`|\chi\rangle\langle\nu\Sigma_f|`
     - A single canonical-polyadic term: one emission column
       :math:`\chi` against one production row :math:`\nu\Sigma_f`. The
       *rank-1 degenerate* of the scattering sum below — fission is the
       :math:`\ell=0`, one-mode case.
   * - **Scattering** :math:`S_{\rm aniso}`
     - orthogonal-CP / spectral sum
       :math:`\sum_\ell |Y_\ell\rangle\,\sigma_\ell\,\langle Y_\ell|`
     - A *sum* of rank-1 dyads in the spherical-harmonic eigenbasis
       (Funk–Hecke): the modes :math:`|Y_\ell\rangle` are orthogonal and
       the weights :math:`\sigma_\ell = \Sigma_{s,\ell}` are the Legendre
       moments — the spectral theorem :math:`U\Sigma U^*`, managed by the
       Frame (:math:`R\circ\Lambda\circ M`).
   * - **Full scattering** :math:`S`
     - block-term decomposition (BTD)
     - rank-1 in **angle** :math:`\otimes` a *full-rank* **energy**
       transfer :math:`\otimes` diagonal in **space**. The energy block
       is genuinely dense (group-to-group transfer is not low-rank);
       fission is its CP-rank-1 degenerate (one energy mode, one angle
       mode).

The collision operator :math:`C = M[\sigma_t]` completes the picture as
the **diagonal** term (:eq:`multiplication-operator-action`, the §5.7
multiplication operator): no decomposition, pointwise in every axis. So
the transport algebra's emission/loss structure reads as one ladder —
diagonal (collision) :math:`\to` rank-1 CP atom (fission) :math:`\to`
orthogonal-CP spectral sum (anisotropic scattering) :math:`\to` block-term
decomposition (full scattering).

.. warning::

   **The decomposition is a lens; the guardrails are load-bearing.** This
   framing must NOT be turned into machinery:

   * **Keep the energy block dense.** The group-to-group transfer is
     genuinely full-rank; do not attempt a low-rank energy factorisation.
   * **Do not import general-CP / tensor-fitting.** ALS / CP-decomposition
     *fitting* (the numerical recovery of unknown factors) is foreign to
     this algebra — the factors here are *known by physics* (:math:`\chi`,
     :math:`\nu\Sigma_f`, the :math:`Y_\ell` eigenbasis), not fitted.
   * **Do not mint a ``CPOperator`` or a polyadic umbrella type.** The dyad
     :func:`~orpheus.numerics.operator.outer`, the
     :class:`~orpheus.transport.frames.harmonic_frame.HarmonicFrame`, and
     ``⊗`` are already the right-grained primitives; an umbrella type would
     add ceremony and a conversion seam without making any illegal state
     unrepresentable (``coding-elegance`` type-vs-property: a separate type
     earns its existence only with :math:`\ge 2` non-isomorphic realisations
     under a non-identity morphism — the dyad and the Frame already *are*
     those two realisations, named honestly).

.. (vv-status rationale) Conceptual lens organising the two emission
   kernels in tensor-decomposition vocabulary; not a type and not a solver
   claim. The verifiable content is the fission dyad's closed-form k∞ gate
   and the scattering kernel's anisotropic MMS gate, both cited above.
.. vv-status: emission-kernels-btd documented

.. _scattering-carrier-grid:

The completed carrier grid — the four leaves and the three edges
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The :math:`R \circ \Lambda \circ M` kernel above reads cleanly as an
:class:`~orpheus.numerics.operator.OperatorProduct` of three ``np.ndarray``
maps, but each factor crosses between two **typed transport carriers** —
the per-ordinate flux and its harmonic moments, in their flux and
source/sink roles. The Frame campaign's P4 phase closed those carriers
into a complete :math:`(\text{angular} \otimes \text{moment}) \times
(\text{flux} \otimes \text{source})` grid: four leaf types and the three
edges that map between them.

.. math::
   :label: scattering-carrier-grid

   \begin{array}{ccc}
     \texttt{AngularFlux}
       & \xrightarrow{\;\;M\;\;}
       & \texttt{HarmonicMomentFlux} \\[2pt]
     & & \big\downarrow{\scriptstyle\;\Lambda} \\[2pt]
     \texttt{AngularSourceSink}
       & \xleftarrow{\;\;R\;\;}
       & \texttt{HarmonicMomentSourceSink}
   \end{array}

.. (vv-status rationale) The carrier-grid layout: a named-field-typing
   identity (which carrier type sits at each node, and the role/axis
   semantics of each edge). Not a solver claim; the verifiable content is
   the role/class-identity algebra of the four leaves (the foundation
   tests ``tests/sn/primitives/test_typed_source_sinks.py`` ::
   ``TestHarmonicMomentSourceSink`` and ``tests/transport/frames/
   test_harmonic_frame.py``) and the 0-ULP equivalence of the typed and
   ndarray scattering arms (``test_scattering_kernel_crosscheck.py``).
.. vv-status: scattering-carrier-grid documented

The **two vertical axes** of the grid are the representation (angular ↔
moment) and the role (flux ↔ source). The four leaves and three edges:

.. list-table:: The four carriers and the three edges of the scattering grid
   :header-rows: 1
   :widths: 22 30 48

   * - Carrier / edge
     - Type
     - What it is
   * - top-left leaf
     - :class:`~orpheus.transport.fields.angular_flux.AngularFlux`
       (``FluxRole``)
     - The per-ordinate flux :math:`\psi_n(\vec r, g)` — an affine flux
       *state* (a point, not a vector; ``flux + flux`` is gated by the
       :class:`~orpheus.transport.fields._flux_role.FluxRole` torsor).
   * - top-right leaf
     - :class:`~orpheus.transport.fields.harmonic_moment_flux.HarmonicMomentFlux`
       (``FluxRole``)
     - The flux moments :math:`\phi_\ell^m(\vec r, g)` — the same flux
       state in moment space (``(L+1, 2L+1, ng, *spatial)``). A
       ``(FluxRole, MomentField)`` carrier; ``flux units`` =
       :data:`~orpheus.numerics.units.SCALAR_FLUX_UNITS` (a moment is
       angle-integrated, so the :math:`\ell=0` block **is** the scalar
       flux exactly).
   * - bottom-right leaf
     - :class:`~orpheus.transport.source_sinks.harmonic_moment_source_sink.HarmonicMomentSourceSink`
       (bare ``MomentField``)
     - The **scattered in-scatter source** moments — the output of
       :math:`\Lambda`. A *rate density*, so it adds vectorially
       (``source + source`` is CLOSED, no affine gate);
       :data:`~orpheus.numerics.units.SCALAR_RATE_UNITS`. The P4 leaf that
       gave the flux→source role change a *home* (before it, the role
       change leaked to the scattering consumer as a raw ``np.ndarray``).
   * - bottom-left leaf
     - :class:`~orpheus.transport.source_sinks.angular_source_sink.AngularSourceSink`
       (bare)
     - The per-ordinate in-scatter source :math:`Q_n(\vec r, g)` the
       sweep consumes — the bottom of the chain.
   * - :math:`M` (left edge)
     - :meth:`HarmonicFrame.analyse <orpheus.transport.frames.harmonic_frame.HarmonicFrame.analyse>`
     - **Role-preserving, axis-changing.** Projects per-ordinate →
       moment, flux→flux *and* source→source (the analysis face of the
       Galerkin frame, the canonical pure-Galerkin :math:`\Pi`).
   * - :math:`\Lambda` (right edge)
     - :meth:`LegendreMomentScattering.apply <orpheus.transport.operators.scattering.LegendreMomentScattering.apply>`
     - **Role-changing, axis-preserving.** The *sole* role-changing edge:
       the per-:math:`\ell` group-transfer :math:`\Sigma_{s,\ell}` maps a
       :class:`HarmonicMomentFlux` to the
       :class:`HarmonicMomentSourceSink` it emits — flux → source, both in
       moment space.
   * - :math:`R` (bottom edge)
     - :meth:`HarmonicFrame.reconstruct <orpheus.transport.frames.harmonic_frame.HarmonicFrame.reconstruct>`
     - **Role-preserving, axis-changing.** Reconstructs moment →
       per-ordinate, flux→flux *and* source→source (the reconstruction
       face; the addition-theorem :math:`R`, not the W-weighted adjoint of
       :math:`M`).

The kernel of the previous subsection is the composite of these three
edges plus the producer-side normalisation,

.. math::

   S_{\rm aniso} \;=\; \tfrac{1}{W}\,(R \circ \Lambda \circ M)
   \;:\; \texttt{AngularFlux} \longrightarrow \texttt{AngularSourceSink},

a flux at the top-left mapped all the way round the grid to the source
at the bottom-left. The role only ever changes **once**, at :math:`\Lambda`
— that is the physical content of "scattering turns flux into source", now
visible in the type signatures rather than buried inside an ndarray chain.
The frame verbs are role-polymorphic by
:func:`~functools.singledispatchmethod`, so the *same* :math:`M` and
:math:`R` carry both the flux leg (top edge / bottom edge, flux side) and
the source leg used by the windowed in-scatter migration below.

This 2×2 scattering square is one face of a larger structure. The next
three subsections lift it to the full :math:`(\text{Representation} \times
\text{Role})` carrier grid, identify that grid as a **double category**
whose 2-cell IS :math:`\texttt{frame.conjugate}(\Lambda)`, and then derive
the load-bearing architectural consequence: why the flat
multiple-inheritance leaves the grid is built from are not a workaround but
the **unique principled normal form**, and why the genericity the grid
expresses belongs on the *operator* (:math:`[\textsf{Domain},
\textsf{Codomain}]`), never on the carrier.


.. _carrier-grid-double-category:

The carrier grid is a double category — two kinds of morphism, one 2-cell
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The scattering square above is not a special diagram — it is one cell of a
grid that every transport carrier and every transport operator inhabits.
A carrier is a pair

.. math::
   :label: carrier-grid-cell

   \texttt{Carrier} \;=\; (\,\text{Representation},\ \text{Role}\,),

and the two coordinates are **independent and orthogonal**, each governing
a different facet of the object:

* **Representation** :math:`\in \{\text{Angular},\ \text{Moment},\
  \text{Scalar},\ \text{Trace}\}` sets the **array shape** and carries the
  change-of-basis. Angular is the per-ordinate :math:`(N, n_g, *\text{spatial})`
  layout; Moment is the harmonic-coefficient
  :math:`(L{+}1, 2L{+}1, n_g, *\text{spatial})` layout; Scalar is the
  angle-integrated :math:`(n_g, *\text{spatial})` layout; Trace is the
  boundary face-cochain (the flat :math:`(\,\text{layout.total\_size}\,)`
  buffer on the :class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace`).
  Changing representation is a change of basis between two *realisations of
  the same physical quantity* — the addition theorem :math:`M`/:math:`R`
  between per-ordinate and moment angular space, the angular integral
  between Angular and Scalar. This axis is the storage-family ABC layer of
  :mod:`orpheus.transport.fields._bases`
  (:class:`~orpheus.transport.fields._bases.AngularField` /
  :class:`~orpheus.transport.fields._bases.MomentField` /
  :class:`~orpheus.transport.fields._bases.ScalarField` /
  :class:`~orpheus.transport.fields._bases.AngularBoundaryField`).

* **Role** :math:`\in \{\text{Flux},\ \text{Source},\ \text{Residual},\
  \text{Displacement}\}` sets the **arithmetic interface** — the #208
  affine torsor of :ref:`affine-typed-field-algebra`. Flux is an affine
  *point* (``flux − flux → Displacement``, ``flux ⊕ Displacement → flux``,
  and ``flux + flux`` is a :class:`TypeError` — :eq:`affine-torsor-algebra`);
  Source, Residual, and Displacement are *vectors* in the associated
  difference space (``source + source`` is closed, ``residual + residual``
  is closed). The role is what the object's ``__add__`` *means*. This axis
  is the role-mixin layer
  (:class:`~orpheus.transport.fields._flux_role.FluxRole` for the flux
  states; the bare storage base, i.e. *no* role mixin, for the Source and
  Residual leaves; :class:`~orpheus.transport.displacements._displacement.Displacement`
  for the increment leaves).

The two axes are genuinely orthogonal: a representation change preserves
role (the addition theorem maps a flux to a flux and a source to a source
— it never turns a flux into a source), and a role change preserves
representation (scattering a moment flux to a moment source stays in
moment space). That orthogonality is exactly the structure of a **double
category**, and naming it that way is not decoration — it tells you which
generic each morphism is, and it identifies a coherence theorem the code
already pins to 0 ULP.

.. list-table:: The carrier grid as a double category
   :header-rows: 1
   :widths: 24 30 46

   * - Categorical part
     - In the carrier grid
     - Consequence for the code
   * - **Objects** (0-cells)
     - The grid cells :math:`(\text{Representation}, \text{Role})` — the
       :math:`\approx 4 \times 4` leaf types
       (:class:`~orpheus.transport.fields.angular_flux.AngularFlux`,
       :class:`~orpheus.transport.source_sinks.harmonic_moment_source_sink.HarmonicMomentSourceSink`,
       …), each a 2-line MI binding ``Leaf(RoleMixin, RepBase)``.
     - There is **no cell-by-cell duplication** — each leaf is the
       intersection of one role mixin and one representation base. The
       grid is already its own normal form (see
       :ref:`carrier-grid-flat-leaf-normal-form`).
   * - **Horizontal 1-morphisms**
     - **Representation-changes** — the frame faces :math:`M` (analysis)
       and :math:`R` (reconstruction), built from the
       :math:`(\text{basis}, \text{measure})` pair that *is* the Frame
       (:ref:`scattering-carrier-grid`). A horizontal arrow fixes the
       Role coordinate.
     - A base change that **fixes the fiber** ⟹ :math:`M`/:math:`R` are
       **role-generic**: the *same* analysis face projects a flux to a
       flux and a source to a source. This is why the frame verbs are
       role-polymorphic, and why the role-changing edge is deliberately
       **not** a frame verb.
   * - **Vertical 1-morphisms**
     - **Role-changes** — the cross sections :math:`C = \sigma_t`
       (collision), :math:`\Lambda = \Sigma_{s,\ell}` (the per-:math:`\ell`
       group transfer), :math:`F = \chi \otimes \nu\Sigma_f` (fission). A
       vertical arrow fixes the Representation coordinate.
     - A fiber morphism **identical over every base** ⟹ the cross sections
       are **representation-generic**: the role change *is* the
       cross-section physics (flux → emitted source), and it carries the
       same meaning whether applied in angular, moment, or scalar
       representation. "Scattering turns flux into source" is a vertical
       arrow.
   * - **The 2-cell**
     - **Scattering** :math:`S_{\rm aniso} = \tfrac{1}{W}\,(R \circ
       \Lambda \circ M)` — the vertical 1-morphism :math:`\Lambda`
       **conjugated by the horizontal adjoint pair** :math:`M`/:math:`R`.
       Realised as :meth:`frame.conjugate(Λ) <orpheus.numerics.frame.FrameBase.conjugate>`
       :math:`=` ``OperatorProduct(R, OperatorProduct(Λ, M))``.
     - The 2-cell fills the square: it is the canonical conjugation of a
       vertical morphism by the horizontal frame. The role change stays
       **localized at** :math:`\Lambda` (:math:`M`/:math:`R` preserve
       role), so :math:`S` is honestly an operator from the
       :math:`(\text{Angular}, \text{Flux})` cell to the
       :math:`(\text{Angular}, \text{Source})` cell.

.. (vv-status rationale) The double-category reading of the carrier grid:
   a structural / categorical identity naming which morphism class each
   operator belongs to (horizontal = representation-change, vertical =
   role-change, scattering = the 2-cell). Not a solver claim; the
   verifiable content is the 0-ULP interchange-coherence identity
   (``tests/sn/operators/test_frame_conjugate_carve.py`` ::
   ``TestFrameConjugateEqualsRLambdaM`` +
   ``test_kernel_property_is_frame_conjugate_of_lambda``) plus the
   role/class-identity algebra of the leaves (the foundation tests cited
   at :eq:`scattering-carrier-grid`).
.. vv-status: carrier-grid-cell documented

The interchange law is a theorem the code already pins
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A double category's **interchange law** states that the 2-cell may be
read either way round the square and the two readings agree. For the
scattering 2-cell this is the statement that the *single composed*
operator :math:`R \circ \Lambda \circ M` and the *step-by-step* typed
evaluation (project with :math:`M`, scatter with :math:`\Lambda`,
reconstruct with :math:`R`) compute the **same** per-ordinate source —
not approximately, but to the last bit. That is exactly what
``tests/sn/operators/test_scattering_kernel_crosscheck.py`` and
``tests/sn/operators/test_frame_conjugate_carve.py`` assert with
``np.array_equal`` (0 ULP, **not** ``allclose``):

.. math::
   :label: carrier-grid-interchange-witness

   \underbrace{\big(\texttt{frame.conjugate}(\Lambda)\big).\texttt{apply}(\psi)}_{\text{single composed 2-cell}}
   \;\;\equiv\;\;
   \underbrace{R\big(\Lambda\,(M\,\psi)\big)}_{\text{step-by-step horizontal·vertical·horizontal}}
   \qquad(0\ \text{ULP}).

.. (vv-status rationale) The interchange-coherence identity: the composed
   2-cell equals the step-by-step horizontal/vertical/horizontal reading
   bit-for-bit. The verifiable content is the 0-ULP ``np.array_equal``
   crosscheck (``test_scattering_kernel_crosscheck.py``, the definitional
   identity of :attr:`ScatteringOperator.kernel`, and
   ``test_frame_conjugate_carve.py`` ::
   ``test_conjugate_equals_manual_R_A_M_nesting``). Not a solver claim —
   a structural equivalence between two evaluations of one operator.
.. vv-status: carrier-grid-interchange-witness documented

The bit-identity is the point. The two sides of
:eq:`carrier-grid-interchange-witness` share the *same* :math:`\Lambda`
kernel and the *same* frame :math:`R` face, so their agreement is a
**coherence theorem of the double category** — it holds by construction,
not by numerical coincidence — and the 0-ULP gate is its
**interchange-law coherence witness**. This is why the crosscheck is a
``np.array_equal`` *definitional* identity rather than an ``allclose``
*regression* tolerance: a tolerance would admit two genuinely different
reduction trees agreeing only to round-off, which would mean the square
does not actually commute. The 0 ULP says it commutes exactly — the mark
of a real 2-cell.

.. note::

   **An equivalent reading: a category fibered over Representation.** The
   same structure can be stated as a *fibered category* (a Grothendieck
   fibration) :math:`p : E \to B` with base :math:`B =` Representation and
   the **Role as the fiber coordinate**, carrying a **torsor on the Flux
   fiber** (:ref:`affine-typed-field-algebra`). In this reading a
   role-change (:math:`\Lambda`, :math:`C`, :math:`F`) is a *cartesian
   morphism within a fiber* (fixed representation), and a
   representation-change (:math:`M`, :math:`R`) is a *base change* lifting
   to the total space. :math:`M` is role-generic **because a base change
   fixes the fiber coordinate and lifts uniformly**; :math:`\Lambda` is
   representation-generic **because the fiber morphism is the same over
   every base point**. The double-category and fibration pictures describe
   the same code; the double category makes the 2-cell (scattering)
   explicit, the fibration makes the role-preservation a *theorem* (base
   change fixes the fiber) rather than a per-operator assertion.


.. _carrier-grid-domain-codomain-identity:

The key identity — a grid cell IS an operator's ``(Domain, Codomain)``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The grid has two coordinates on the *carrier*; an operator has two
coordinates too — its **domain** and its **codomain**. These are the same
two coordinates. An operator is an arrow from one grid cell to another,
and naming its endpoints names the morphism completely:

.. math::
   :label: carrier-grid-operator-typing

   \texttt{LinearOperator[Domain, Codomain]}
   \;\;\text{IS the typed traversal of the grid:}\quad
   \begin{aligned}
     M &: \texttt{LinearOperator[AngularFlux,\ HarmonicMomentFlux]}
        &&\text{(horizontal)}\\
     \Lambda &: \texttt{LinearOperator[HarmonicMomentFlux,\ HarmonicMomentSourceSink]}
        &&\text{(vertical)}\\
     S &: \texttt{LinearOperator[AngularFlux,\ AngularSourceSink]}
        &&\text{(the 2-cell)}.
   \end{aligned}

.. (vv-status rationale) The operator-typing identity: an operator's two
   type parameters ARE the two grid cells it maps between. A
   representational/structural statement about where the parametrization
   lives (on the operator, not the carrier). The verifiable content is the
   static ``assert_type`` pins on the heteromorphic ``apply``
   (``tests/sn/operators/test_operators_apply_typed.py``) and the
   composition-guard tests; not a solver claim.
.. vv-status: carrier-grid-operator-typing documented

The two-parameter operator type
:class:`~orpheus.numerics.operator.LinearOperator` (``Protocol[Domain,
Codomain]``, :ref:`heteromorphic-apply-typing`) is therefore the **right
and complete machinery for traversing the grid**. Its ``apply`` maps an
input carrier :data:`~orpheus.numerics.operator.Domain` to a
(possibly distinct) output carrier
:data:`~orpheus.numerics.operator.Codomain` — exactly an arrow
:math:`(\text{Rep}_{\rm in}, \text{Role}_{\rm in}) \to (\text{Rep}_{\rm
out}, \text{Role}_{\rm out})`. The endomorphic majority (collision
:math:`C`, identity, the ``np.ndarray`` serialization boundary) is the
special case :math:`\textsf{Codomain} = \textsf{Domain}`, recovered for
free by the PEP-696 default ``Codomain = TypeVar("Codomain",
default=Domain)`` so that ``LinearOperator[V]`` :math:`\equiv`
``LinearOperator[V, V]`` and the endomorphic call sites need no change.

The names are spelled in full — :data:`~orpheus.numerics.operator.Domain`
and :data:`~orpheus.numerics.operator.Codomain`, not abbreviations —
because ``Domain`` already reads as "the input" and ``Codomain`` as "the
output"; the morphism vocabulary is the domain vocabulary here.

.. important::

   **The parametrization belongs on the operator, NOT on the carrier.**
   This is the load-bearing architectural decision the grid forces. One
   could imagine pushing the :math:`(\text{Representation}, \text{Role})`
   pair *onto the carrier* as type parameters —
   ``Carrier[Representation, Role]`` — and writing operators generic in
   both axes. That design is **structurally impossible in Python** and
   would break a runtime safety gate even where it type-checks; the next
   subsection is the full argument. The realised design puts the two
   coordinates where they are expressible and where they belong: on the
   **operator's** ``[Domain, Codomain]``, with the carriers as the flat
   intersection leaves the grid cells already are.


The HarmonicFrame typed seam — and why it lives in transport, not numerics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The casting between the generic ``np.ndarray`` frame faces and the typed
:class:`Field` carriers is the load-bearing design decision of P4, and it
is forced to live **one layer up from the frame itself**. The argument is
a layering constraint, not a preference:

* The generic angular spherical-harmonic frame is the numerics
  :class:`~orpheus.numerics.frame.GalerkinFrame` built by
  :meth:`~orpheus.numerics.quadrature.Quadrature.angular_frame` — its
  :attr:`~orpheus.numerics.frame.FrameBase.analysis` /
  :attr:`~orpheus.numerics.frame.FrameBase.reconstruction` faces are
  **carrier-agnostic** ``np.ndarray → np.ndarray`` maps. This is by
  design: the same generic faces are shared with the P3
  indicator-homogenisation frame, and keeping them ndarray-valued is what
  makes that sharing 0-ULP-safe.
* The two carriers the angular frame maps between
  (:class:`AngularFlux` ↔ :class:`HarmonicMomentFlux`, and their
  source/sink siblings) share their *deepest* primitive,
  :class:`~orpheus.numerics.field.Field`, in **numerics**.
* But the part that makes them **castable** — the ``mesh`` binding plus
  the ``from_mesh`` / ``from_mesh_and_L`` factories that build the typed
  carrier from a raw array — lives in the transport
  :class:`~orpheus.transport.fields._bases.BulkField` base, **above**
  numerics. And :meth:`Quadrature.angular_frame <orpheus.numerics.quadrature.Quadrature.angular_frame>`
  is in numerics, which *cannot* import the transport carriers without
  inverting the layer order.

So a generic numerics face **cannot** return a typed
:class:`HarmonicMomentFlux`: the cast needs the transport-layer factories,
and numerics is below transport. The clean home for the casting is
therefore the transport layer, in
:class:`~orpheus.transport.frames.harmonic_frame.HarmonicFrame`:

.. math::

   \texttt{HarmonicFrame} \;\text{IS-A}\; \texttt{GalerkinFrame}
   \qquad(\text{Liskov — the angular SH projection is the canonical
   pure-Galerkin frame}),

constructed from the generic frame's basis + measure via
:meth:`~orpheus.transport.frames.harmonic_frame.HarmonicFrame.from_galerkin`
(``cls(frame.basis, frame.measure)`` — **zero** rebuild of the basis /
measure / projection table, so the inherited ndarray faces and the §5.6
:attr:`kernel <orpheus.transport.operators.scattering.ScatteringOperator.kernel>` stay
bit-identical), adding **only** the typed verbs
:meth:`~orpheus.transport.frames.harmonic_frame.HarmonicFrame.analyse`
(:math:`M`) and
:meth:`~orpheus.transport.frames.harmonic_frame.HarmonicFrame.reconstruct`
(:math:`R`). The generic numerics faces are untouched. This is the
"shared primitive is :class:`Field` in numerics, but castability is
:class:`BulkField` in transport, so the typed seam lives in transport"
layering — the same reason a typed wrapper, not the generic frame, owns
the carrier verbs.

.. note::

   The role-changing edge :math:`\Lambda` is deliberately **not** a frame
   verb. A frame's faces are role-*preserving* changes of representation
   (flux↔flux, source↔source); the flux→source role change is *physics*
   (scattering emits a source), so it lives on the scattering operator
   (:meth:`LegendreMomentScattering.apply <orpheus.transport.operators.scattering.LegendreMomentScattering.apply>`),
   where the cross sections are. Putting :math:`\Lambda` on the frame would
   conflate "change of basis" with "change of physical kind".

Explicit typed path vs the fused composed kernel — the windowed arm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

P4 deliberately keeps **two** realisations of the same :math:`R\Lambda M`
math, each chosen for what its consumer needs to see:

* The hot **full-angular** path (:ref:`sn-angular-windowing-factoring`,
  :meth:`~orpheus.transport.operators.scattering.ScatteringOperator.build_aniso_source`)
  consumes the §5.6 :attr:`kernel <orpheus.transport.operators.scattering.ScatteringOperator.kernel>`
  ``= frame.conjugate(Λ)`` — the **single composed** ``np.ndarray``
  operator. This is P2's win and the 0-ULP canary: one composition, one
  reduction tree, pinned bit-for-bit by
  ``tests/sn/operators/test_scattering_kernel_crosscheck.py``. The role
  change is *implicit* inside the ndarray chain; that is correct here,
  because the consumer is a tight numerical loop that never names the
  intermediate moment source.
* The **windowed** in-scatter arm (the
  :class:`~orpheus.transport.fields.harmonic_moment_flux.HarmonicMomentFlux`
  arm of
  :meth:`ScatteringOperator.apply <orpheus.transport.operators.scattering.ScatteringOperator.apply>`)
  takes the **explicit typed** edges: :math:`\Lambda` scatters the flux
  moments to a typed :class:`HarmonicMomentSourceSink` (the role-changing
  edge made visible in the signature), then
  :meth:`frame.reconstruct <orpheus.transport.frames.harmonic_frame.HarmonicFrame.reconstruct>`
  reconstructs the per-ordinate :class:`AngularSourceSink`, then the
  producer-side :math:`1/W`. The realisation is
  ``frame.reconstruct(Λ.apply(φ)) / W`` — the explicit, role-typed
  counterpart of the fused kernel.

Both arms route through the *same* :math:`\Lambda` kernel and the *same*
frame :math:`R` face, so they agree numerically. The choice is one of
**legibility at the call site**, not of math: the windowed consumer holds
the moments as a typed iterate and the explicit edges make the
flux→source→angular role flow read off the signatures; the hot consumer
holds a raw ndarray and the fused operator keeps the reduction tree single
and the 0-ULP canary meaningful. The ndarray ``reconstruct_after(Λ)``
reference
(:meth:`~orpheus.transport.operators.scattering.ScatteringOperator._aniso_source_from_moment_values`)
is retained as the crosscheck's oracle — the structurally-independent
``np.ndarray`` evaluation the typed arm is pinned against.


.. _carrier-grid-census:

The full (Representation × Role) census — and its two principled holes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The scattering square uses four cells of the grid; the whole grid is the
:math:`4 \times 4` product of the four representations against the four
roles. Most cells are realised as a concrete leaf; two are *deliberately*
empty, for reasons that are themselves part of the architecture (an empty
cell that is principled tells a future session "do not mint this", which is
as load-bearing as a populated one). The census below is the live state of
:mod:`orpheus.transport.fields`,
:mod:`orpheus.transport.source_sinks`,
:mod:`orpheus.transport.residuals`, and
:mod:`orpheus.transport.displacements`.

.. list-table:: The (Representation × Role) carrier grid — leaf census
   :header-rows: 1
   :stub-columns: 1
   :widths: 18 21 21 21 19

   * -
     - **Flux** (``FluxRole``)
     - **Source** (bare)
     - **Residual** (bare)
     - **Displacement**
   * - **Angular**
     - :class:`~orpheus.transport.fields.angular_flux.AngularFlux`
     - :class:`~orpheus.transport.source_sinks.angular_source_sink.AngularSourceSink`
     - :class:`~orpheus.transport.residuals.angular_residual.AngularResidual`
     - :class:`~orpheus.transport.displacements.angular_displacement.AngularDisplacement`
   * - **Moment**
     - :class:`~orpheus.transport.fields.harmonic_moment_flux.HarmonicMomentFlux`
     - :class:`~orpheus.transport.source_sinks.harmonic_moment_source_sink.HarmonicMomentSourceSink`
     - — *(principled hole, below)*
     - :class:`~orpheus.transport.displacements.moment_displacement.MomentDisplacement`
   * - **Scalar**
     - :class:`~orpheus.transport.fields.scalar_flux.ScalarFlux`
     - :class:`~orpheus.transport.source_sinks.scalar_source_sink.ScalarSourceSink`
     - :class:`~orpheus.transport.residuals.scalar_residual.ScalarResidual`
     - :class:`~orpheus.transport.displacements.scalar_displacement.ScalarDisplacement`
   * - **Trace**
     - :class:`~orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux`
     - :class:`~orpheus.transport.source_sinks.angular_boundary_source_sink.AngularBoundarySourceSink`
     - :class:`~orpheus.transport.residuals.angular_boundary_residual.AngularBoundaryResidual`
     - :class:`~orpheus.transport.displacements.angular_boundary_displacement.AngularBoundaryDisplacement`

Reading the columns confirms the two-axis structure of
:ref:`affine-typed-field-algebra`: the **Flux** column carries the
:class:`~orpheus.transport.fields._flux_role.FluxRole` torsor mixin (an
affine point); the **Source** and **Residual** columns are *bare* storage
leaves (plain vector algebra — ``source + source`` and ``residual +
residual`` are closed, no affine gate); the **Displacement** column carries
the :class:`~orpheus.transport.displacements._displacement.Displacement`
mixin (the iterate-increment vector, home of the contraction diagnostics).
The Trace (boundary) row mirrors the bulk rows exactly — the four-column
parallel the boundary role grid completes at
:ref:`bc-extraction-operator-output-typing`.

The role-axis "asymmetry" is the type-vs-property rule, not a defect
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A first glance at the columns sees an asymmetry: Flux and Displacement are
**mixins** (they add behaviour — the torsor algebra, the contraction
diagnostics), while Source and Residual carry **no mixin** (they are bare
representation leaves). This is **correct**, and it is the project's
type-vs-property rule (CLAUDE.md Cardinal Rule 2; the
``coding-elegance`` "build the primitive" pattern) applied exactly: mint a
distinct **role object** only where a **non-identity morphism** lives on
that role.

* **Flux** earns a class because it has special arithmetic — the torsor
  (``flux + flux`` must *raise*, ``flux − flux`` must mint a
  *Displacement*). That behaviour is unrepresentable without a class to
  carry the overridden ``__add__`` / ``__sub__``.

* **Displacement** earns a class because it carries diagnostics a flux
  state structurally cannot (``contraction_ratio``, the Aitken / true-error
  estimate — :ref:`affine-typed-field-algebra`).

* **Source** and **Residual** are *plain vector roles* — their arithmetic
  is the inherited :class:`~orpheus.numerics.field.Field` vector algebra,
  with **no** special dunder. There is no non-identity morphism to host, so
  by the rule they **correctly carry nothing** beyond their representation
  base. The role distinction between a Source and a Residual is still real
  — they are different *classes* (a residual is born only from a
  :meth:`~orpheus.numerics.field.Field._from_balance`, a source from an
  operator ``apply``; the class identity gates cross-role addition even
  though they share rate-density units), but the distinction needs no
  *behaviour*, so neither gets a mixin.

Uniformising the role axis — giving Source and Residual an empty marker
mixin "for symmetry" — would be **ceremony**: it would add a class with no
behaviour, exactly the kind of theatrics the type-vs-property rule forbids.
The asymmetry is the rule working.

The two principled holes
^^^^^^^^^^^^^^^^^^^^^^^^^

Two cells are deliberately empty, and both absences are designed.

* **(Moment, Residual) —** :class:`!HarmonicMomentResidual` **is absent.**
  A residual is born only from a balance equation
  (:meth:`~orpheus.numerics.field.Field._from_balance`), and **moment space
  is never the subject of a balance**. The transport balances are the
  bulk-angular :math:`(L + C - S - F)\psi - q` and the boundary consistency
  :math:`\psi.\text{inflow} - B\,\psi.\text{outflow} - q.\text{inflow}`
  (:ref:`bc-extraction`); neither is posed in moment coordinates, so there
  is no ``from_balance`` consumer that would produce a moment residual.
  Minting :class:`!HarmonicMomentResidual` would create a leaf no producer
  fills — an illegal state by the "build the primitive only when it has a
  consumer" rule. The hole is the rule, not an oversight. (Contrast
  **(Moment, Displacement)**, which *is* populated:
  :class:`~orpheus.transport.displacements.moment_displacement.MomentDisplacement`
  exists because the angular-windowed SI iterate *does* hold its state as a
  moment flux and *does* difference it — :ref:`sn-angular-windowing` — so
  the moment displacement has a genuine consumer.)

* **The** ``iso + aniso → AngularSourceSink`` **source injection is a
  hand-rolled Representation traversal inside a dunder, and it is
  endorsed.** :meth:`AngularSourceSink.__add__ <orpheus.transport.source_sinks.angular_source_sink.AngularSourceSink.__add__>`
  accepts a :class:`~orpheus.transport.source_sinks.scalar_source_sink.ScalarSourceSink`
  partner and applies the **canonical subspace-containment injection**: a
  scalar (isotropic) source lives in the subspace of the per-ordinate
  (angular) source where every ordinate carries the same value, and the
  injection :math:`\text{iso} \to \mathbf 1 \otimes \text{iso}` (broadcast
  across the :math:`\Omega` axis) maps it in before the add. This is a
  *Representation* change (Scalar → Angular) performed inside a *Role*
  operation (source + source), so it sits slightly outside the clean
  "horizontal morphisms are frame faces" story — it is a one-off
  representation embedding baked into a leaf's arithmetic. It is kept
  because the embedding is the genuine mathematical relation between an
  isotropic and an anisotropic source (the :math:`\ell=0` block *is* the
  isotropic part), it is single-sourced through
  :meth:`~orpheus.transport.source_sinks.angular_source_sink.AngularSourceSink.from_isotropic`,
  and it was **endorsed at creation** as the right home for the iso/aniso
  combine. It is documented here as a recognised, deliberate exception so a
  future reader does not "tidy it away".


.. _carrier-grid-flat-leaf-normal-form:

Why a fully-typed ``Carrier[Representation, Role]`` is impossible
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The grid *invites* an obvious design: make the carrier itself generic in
both axes — ``Carrier[Representation, Role]`` — and write operators that
are generic in both, so that the whole grid is one parametrized type and
the leaf zoo collapses. This is the design the architecture **explored and
rejected**, because it is structurally impossible in Python and would break
a runtime safety gate even in the fragments where it type-checks. Recording
*why* it cannot work is the point of this subsection: it stops a future
session from re-attempting the "obvious" collapse and discovering the wall
the hard way. The flat multiple-inheritance leaves are not a compromise
forced by a weak type system — they are the **unique principled normal
form**, and the genericity the grid expresses has a correct home (the
operator's ``[Domain, Codomain]``) that the carrier does not.

The argument is five obstructions, each fatal on its own.

.. list-table:: Why ``Carrier[Representation, Role]`` cannot be built
   :header-rows: 1
   :widths: 8 30 62

   * - #
     - Obstruction
     - Why it is fatal
   * - **(a)**
     - **Role changes the arithmetic interface ⟹ Role MUST be a class.**
     - A phantom type parameter (``Generic[Role]``) is **erased at
       runtime**, so it cannot specialize a dunder. But the role *is* the
       ``__add__`` semantics: the Flux role must make ``flux + flux``
       **raise** (the torsor, :eq:`affine-torsor-algebra`) while the Source
       role must make ``source + source`` **succeed**. A single
       ``__add__`` body shared across ``Carrier[Rep, Flux]`` and
       ``Carrier[Rep, Source]`` under one phantom ``Role`` **cannot** make
       one raise and the other succeed — the parameter that would select
       between them does not exist at the moment ``__add__`` runs. This
       hard-refutes the phantom-``Role`` encodings (``Field[Rep, Role]``
       and the representation-outer ``Angular[Role]``). A *runtime*
       ``role`` field that ``__add__`` branches on is the stringly-typed
       anti-pattern (an illegal state is representable —
       ``replace(f, role=Source)`` would bypass the torsor gate), so that
       escape is closed too.
   * - **(b)**
     - **Representation changes the array shape ⟹ Representation MUST be a
       class.**
     - The representation sets the ``values`` / ``space`` shape
       (:math:`(N, n_g, *\text{spatial})` vs
       :math:`(L{+}1, 2L{+}1, n_g, *\text{spatial})` vs the flat trace
       buffer) and the shape-validation hook
       (:meth:`~orpheus.transport.fields._bases.BulkField._phase_space_shape`).
       A phantom ``Representation`` parameter, erased, carries none of that
       — the shape check has nothing to read. This refutes the role-outer
       phantom ``Flux[Rep]`` (the one encoding the torsor obstacle would
       have spared — Role is the outer class there — but the shape
       obstacle kills it instead).
   * - **(c)**
     - **The only both-classes form with role-arithmetic-once and
       rep-shape-once is the flat MI leaf — which already exists.**
     - (a) forces Role to be a class; (b) forces Representation to be a
       class. The form that has *both* as classes, with the role algebra
       written once (per role mixin) and the representation shape written
       once (per storage base), and no per-cell duplication, is the
       multiple-inheritance intersection ``AngularFlux(FluxRole,
       AngularField)`` the grid is **already** built from. There is no
       novel encoding to discover; the normal form is the current code.
   * - **(d)**
     - **A parameterized carrier would break the runtime units gate via
       erasure.**
     - The field algebra enforces units/meaning at runtime with a
       **class-identity** check — ``type(self) is type(other)`` in
       :meth:`~orpheus.numerics.field.Field._check_partner` (extended by
       :meth:`BulkField._check_partner <orpheus.transport.fields._bases.BulkField._check_partner>`
       for the mesh binding). Under a generic ``Carrier[Rep, Role]`` the
       runtime class is the *erased* ``Carrier`` for **every** cell, so the
       identity check would read all cells as the same type and **admit**
       cross-representation, cross-role addition that is physically
       meaningless (adding a moment flux to an angular source). The generic
       carrier does not merely fail to *help* — it actively **disables** a
       working safety gate. The flat leaves keep one concrete class per
       cell, so the identity check stays sharp.
   * - **(e)**
     - **Both-axes-generic operators would need higher-kinded types, which
       Python lacks.**
     - Even granting a generic carrier, an operator generic in *both* axes
       — "for any Representation :math:`X` and any Role :math:`Y`,
       :math:`\Lambda` maps ``Carrier[X, Flux]`` to ``Carrier[X,
       Source]``" — quantifies over a *type constructor* applied to a
       parameter, i.e. a higher-kinded type. Python's type system has no
       higher-kinded type variables, so this signature is unspellable.

**Conclusion.** The flat multiple-inheritance leaves are the **unique,
principled normal form** — not a workaround the language imposes, but the
one design that (i) keeps the torsor's ``flux + flux`` ban expressible,
(ii) keeps the per-representation shape check expressible, (iii) keeps the
runtime units gate sharp, and (iv) avoids per-cell duplication. The
genericity the grid genuinely has — "an operator maps one cell to another"
— lives where Python *can* express it and where it belongs: on the
**operator's** :math:`[\textsf{Domain}, \textsf{Codomain}]`
(:eq:`carrier-grid-operator-typing`), as the typed traversal of the grid.
The carrier stays a flat leaf; the operator carries the two coordinates.

.. note::

   **This is a closed exploration, recorded so it is not re-opened.** The
   four candidate carrier encodings — a phantom ``Field[Representation,
   Role]``, a role-outer ``Flux[Representation]``, a representation-outer
   ``Angular[Role]``, and the current flat MI leaves — were weighed against
   the obstructions above. Only the flat MI leaves survive: the
   phantom-``Role`` forms die on (a), the role-outer form dies on (b), and
   any parameterized carrier dies on (d). The full structural verdict (the
   double-category / fibration / torsor frame attack that produced this
   conclusion) is recorded in
   ``.claude/agent-memory/cross-domain-attacker/rep_role_grid_double_category_frames.md``.
   A future session reaching for ``Carrier[Rep, Role]`` should read this
   subsection first.


Deferred relocation
-------------------

The fission and scattering operators are reframed *in place* in ``sn/``;
their carrier-agnostic, cross-section-only **cores** (the rank-1
production core, the :math:`\Lambda` group-transfer, the iso / :math:`(n,2n)`
fast paths) are the natural residents of the shared ``transport/`` layer
(L2), with only the quadrature-coupled per-ordinate angular adapters
staying in ``sn/`` (L3). That relocation — together with the CP / MoC
carrier unification that would let those solvers *consume* the shared
cores instead of reimplementing C / F / S inline on bare scalar arrays —
is tracked on #261. The dispatch-spelling decision for the heteromorphic
``apply`` (Pattern M vs an ``@overload`` + ``match`` router over the
relocated primitives, see :ref:`heteromorphic-apply-typing`) is parked on
the same issue, to be settled *with* the relocation, because the
router-over-shared-primitives shape falls out of the sharing.

.. note::

   **Source map.** Category Protocol:
   :class:`orpheus.transport.operators.integral_kernel_operator.IntegralKernelOperator`
   (L2). Named kernels:
   :attr:`orpheus.transport.operators.fission.FissionOperator.kernel` +
   :attr:`orpheus.transport.operators.fission.FissionOperator.production_rate`;
   :attr:`orpheus.transport.operators.scattering.ScatteringOperator.kernel`.
   Category-refinement gate:
   ``tests/transport/test_integral_kernel_category.py``. Fission
   cross-check: ``tests/sn/operators/test_fission_kernel_crosscheck.py``
   (hand-derived correctness reference + the
   :math:`\chi \cdot \mathrm{production\_rate} \equiv F.\mathrm{apply}`
   0-ULP de-risk). Scattering cross-check:
   ``tests/sn/operators/test_scattering_kernel_crosscheck.py`` (the
   :math:`S.\mathrm{kernel}.\mathrm{apply} \equiv R\circ\Lambda\circ M`
   0-ULP equivalence). Deferred follow-ups: #260
   (:class:`~orpheus.numerics.operator.SumOfTensorProductsOperator`
   un-orphaning), #261 (core relocation + CP / MoC carrier unification).


.. _tensorial-framing:

Tensor product algebra
======================

Streaming, scattering, and any operator on a multi-axis flux
factor naturally as a **tensor product** of per-axis operators:

* **Streaming** (Grand Report v3 §15.1, line 2044):

  .. math::
     :label: streaming-as-tensor-product-sum

     L \;=\; D_x \otimes \Omega_x \otimes I_g
            + D_y \otimes \Omega_y \otimes I_g.

* **Pℓ moment scattering** (§15.2):

  .. math::
     :label: scattering-as-tensor-product-sum

     S \;=\; \sum_{\ell} P_\ell \otimes \Sigma_{s,\ell},

  where :math:`P_\ell` selects the :math:`\ell`-block on the
  harmonic-coefficient axis and :math:`\Sigma_{s,\ell}` is the
  per-:math:`\ell` group-to-group transfer matrix.

.. vv-status: streaming-as-tensor-product-sum documented
.. vv-status: scattering-as-tensor-product-sum documented

These canonical forms motivate the
:class:`~orpheus.numerics.operator.TensorProductOperator` and
:class:`~orpheus.numerics.operator.SumOfTensorProductsOperator`
types. The user's architectural correction that drove their
introduction:

  *"Is there any tensorial machinery to be brought in support of
  this?"*

The answer is yes, and the cost of NOT having it is that
:math:`M`, :math:`R`, :math:`\Lambda` would be three named
operators that *happen* to be expressible via
:func:`numpy.einsum`, with the axis structure implicit inside each
op. With the tensor-product types, the §9 literal statement
:math:`S_{\text{SN}} = R \Lambda M` becomes the operator-algebra
type signature, and the multi-axis structure is **structural**
rather than buried in einsum subscripts.


Definition
----------

For operators :math:`A_1, A_2, \ldots, A_k` acting on
**independent** tensor axes (each carries an ``axis`` attribute and
broadcasts on the rest), the tensor-product operator's action is
the sequential per-axis application

.. (vv-status rationale) Verified by
   ``tests/numerics/test_tensor_product_operator.py`` — Kronecker-
   product reference on small concrete factors, the
   ``is_invertible`` / ``is_adjointable`` predicate meet, and the
   algebraic laws below.
.. vv-status: tensor-product-action documented

.. math::
   :label: tensor-product-action

   (A_1 \otimes A_2 \otimes \cdots \otimes A_k)\,x
   \;=\; A_k\bigl(\cdots A_2(A_1\,x) \cdots\bigr).

Because the constituents act on disjoint axes, the order does not
matter — the operators commute on the joint tensor. The two-axis
predicates are the **meet** (recursive ``and``) over the
constituents: the tensor product is invertible iff every factor is
invertible, adjointable iff every factor is adjointable (``apply`` is
universal).


Algebraic laws
--------------

The
:class:`~orpheus.numerics.operator.TensorProductOperator` carries
three algebraic laws verified by tests
(:file:`tests/numerics/test_tensor_product_operator.py`):

.. math::
   :label: tensor-product-adjoint-distributivity

   (A \otimes B)^*
   \;=\; A^* \otimes B^*.

.. math::
   :label: tensor-product-axis-wise-composition

   (A \otimes B) \circ (C \otimes D)
   \;=\; (A \circ C) \otimes (B \circ D)
   \quad
   \text{when } A, C \text{ share an axis and } B, D \text{ share
   an axis}.

.. math::
   :label: tensor-product-inverse

   (A \otimes B)^{-1}
   \;=\; A^{-1} \otimes B^{-1}
   \quad
   \text{when both factors are invertible}.

.. vv-status: tensor-product-adjoint-distributivity documented
.. vv-status: tensor-product-axis-wise-composition documented
.. vv-status: tensor-product-inverse documented


Relation to numpy primitives
-----------------------------

This is the load-bearing distinction the user explicitly asked
about:
:func:`numpy.kron`, :func:`numpy.tensordot`, and
:func:`numpy.einsum` are **array primitives** — the
*implementation* layer.
:class:`TensorProductOperator` is the **operator-algebra type** —
a layer above.

.. list-table:: Two abstraction layers
   :header-rows: 1
   :widths: 26 36 38

   * - Layer
     - Carrier
     - Examples
   * - **Implementation** (numpy)
     - Untyped multi-axis arrays; subscript strings encode axis
       structure; algebra is in the programmer's head.
     - ``np.einsum("nlm,n,n...->lm...", Y, w, psi)`` — the axes
       ``n``, ``l``, ``m``, and ``...`` are conventions only;
       nothing in the type system prevents passing a wrong-shaped
       ``Y`` or a wrong-axis ``w``.
   * - **Operator algebra** (this module)
     - Typed operators; ``axis`` attribute on each factor; the two-axis
       predicate meet; algebraic laws checked at composition.
     - ``quad.angular_frame(L).analysis & IdentityOperator()`` —
       the type signal carries the axis structure; mismatched axes
       raise at composition; the :meth:`apply` routes through
       ``np.einsum`` internally with the correct subscripts.

The two layers are **complementary**, not competitive. The
operator-algebra layer routes through the array layer for
performance — every :meth:`apply` call eventually calls
``np.einsum`` or a broadcast-multiply, because numpy's einsum
backend is more optimised than anything Python-level the project
could write. The advantage of the operator-algebra layer is in the
**composition language**:

* ``A & B & C`` reads as :math:`A \otimes B \otimes C`.
* ``(A & B) @ (C & D)`` automatically reduces to
  ``(A @ C) & (B @ D)`` (the axis-wise composition law).
* ``(A & B).H`` is ``A.H & B.H`` (adjoint distributivity).
* The ``is_invertible`` / ``is_adjointable`` predicates are computed
  automatically; a composition mismatch raises at composition, not at
  the first :func:`numpy.einsum` call mid-iteration.

.. note::

   :func:`numpy.kron` constructs a **dense** Kronecker product
   matrix — the explicit :math:`(N_A N_B) \times (N_A N_B)` entry
   table. This is the wrong abstraction for transport: an SN
   sweep's effective tensor-product operator
   :math:`L = D_x \otimes \Omega_x \otimes I_g` would be a
   :math:`(N N_x N_y N_g)^2` matrix — never materialised in
   practice because the action is matrix-free. The
   :class:`TensorProductOperator` carries the algebra without
   the materialisation, just as
   :class:`scipy.sparse.linalg.LinearOperator` carries dense-matrix
   algebra without the dense matrix.


SumOfTensorProductsOperator
---------------------------

When several tensor products are summed,
:class:`~orpheus.numerics.operator.SumOfTensorProductsOperator`
provides the §15.2 canonical scattering / streaming form's named
type:

.. math::
   :label: sum-of-tensor-products

   T \;=\; \sum_{k=1}^{K} A_k \otimes B_k \otimes C_k.

.. vv-status: sum-of-tensor-products documented

Algebraically just :class:`OperatorSum` over
:class:`TensorProductOperator` summands, but exposed as a named
class so the §15.2 invariants
(:meth:`assert_separable` — every summand is a tensor product;
shared-axis-factor refactoring; future TT-compression entry point
per §15.3) carry a load-bearing type signal. The
:meth:`assert_separable` method is currently a contract-validator
(separability is enforced at construction); it is a hook for
future invariant checks.


Tensor product as the inverse of partition
-------------------------------------------

The :class:`~orpheus.numerics.measure.DiscreteMeasure`'s
:meth:`partition_by` method (see :ref:`discrete-measures`) realises
the **direct sum**
:math:`\mu_{S^2} = \bigoplus_\lambda \mu_\lambda`. When the
predicate is the octant-sign label
:math:`\lambda(\hat\Omega) =
(\mathrm{sign}\,\mu_x, \mathrm{sign}\,\mu_y, \mathrm{sign}\,\mu_z)`,
the partition recovers the eight octants of :math:`S^2` (or four
in 2-D where :math:`\mu_z = 0` is a degenerate case).

For per-octant operators (e.g. the SN sweep's :math:`A_{oct}^{-1}`
acting on
:math:`(\text{octant\_ordinates} \times \text{cells} \times
\text{groups})`), the **tensor product** factors the per-axis
structure within an octant while the **direct sum** assembles the
octants into the global angular cubature:

.. math::
   :label: octant-direct-sum-tensor-product

   A^{-1} \;=\; \bigoplus_{\text{oct}}\,
                A_{oct}^{-1}, \qquad
   A_{oct}^{-1} \;=\; \text{(per-axis tensor product within an octant)}.

.. vv-status: octant-direct-sum-tensor-product documented

This is the operator-algebra type signature behind Wave 2 of the
SN performance plan: the 2-D Cartesian sweep iterates
:math:`\bigoplus_{oct}` (4 octants — structural) and per-octant
anti-diagonals (sweep-DAG topology — structural), with no Python
loop over the ordinate axis (which is **internal** to every
:meth:`apply` call within an octant, vectorised by the tensor-
product structure).


.. _field-type-vs-property-criterion:

When a moment representation earns a type (#263)
================================================

The tensor-product algebra above raises a recurring design question whenever
a new moment representation appears: should it be a first-class **field type**
(a sibling of
:class:`~orpheus.transport.fields.harmonic_moment_flux.HarmonicMomentFlux`),
or merely a **property** — a moment axis riding on an existing field?  The
question surfaced sharpest in the SN linear-discontinuous (LD) boundary work
(:ref:`ld-cartesian-2d-coherent-promise`, Issue #257 S9), which had to decide
whether the transverse boundary moment deserved a ``BoundaryMomentField`` type.
This section records the criterion that answers it, because the answer is a
durable design invariant, not a one-off call.

The criterion: a non-canonical dual must coexist
------------------------------------------------

   A representation earns a distinct first-class **type** if and only if there
   exist **two bases that are NOT canonically isomorphic** (the isomorphism
   depends on a quadrature or node choice), connected by a **change-of-basis
   operator that is itself modelled and applied** — it carries truncation
   error, has an adjoint, and participates in the operator algebra.

All three clauses must hold.  This is the sharp form of "a dual must coexist
and not mix", and it is decidable by inspection: count the within-axis
representations and count the applied, non-identity morphisms between them.  If
there is one representation, or the only morphism is the identity, the
representation is a **property**; a type would add no behaviour beyond class
identity — type-theatrics by the project's own standard (a type hint that does
not prevent a bug by construction earns nothing).

The criterion is the field-type analogue of the tensor-product algebra's own
discipline (a typed operator is justified when it carries an ``axis`` attribute,
the two-axis inverse/adjoint predicates, and algebraic laws checked at
composition — not merely a name):
a typed field is justified when it carries a *dual* whose mixing must be
forbidden, not merely a name.

Angular order PASSES; spatial order FAILS (today)
-------------------------------------------------

**Angular order is correctly TWO types.** The ordinate basis
(:class:`~orpheus.transport.fields.angular_flux.AngularFlux`, :math:`N`
collocation directions on :math:`S^2`) and the harmonic-modal basis
(:class:`~orpheus.transport.fields.harmonic_moment_flux.HarmonicMomentFlux`,
:math:`(L+1)(2L+1)` real-spherical-harmonic coefficients) are NOT canonically
isomorphic — the isomorphism depends on the quadrature
:math:`Y_\ell^m(\hat\Omega_n)`.  They are bridged by the APPLIED
projection / reconstruction pair :math:`M` / :math:`R`
(:eq:`harmonic-moment-projection` and its reconstruction inverse), which carry
truncation content and have adjoints and live in the operator algebra.  All
three clauses hold, so the two field types are load-bearing: a ``flux +
moments`` addition is type-rejected by construction (the field-layer partner
gate), exactly as it should be — the ordinate and modal representations must
not silently mix.

**Spatial order is correctly a PROPERTY (today).** There is ONE within-cell
spatial basis — the tensor-Legendre DG tower
(:class:`~orpheus.numerics.spaces.spatial_moment_space.SpatialMomentSpace`,
``per_axis**ndim`` coefficients).  The only change-of-basis within it is the
identity (and ``truncate`` / inclusion, which stay within the same family and
return the same tower).  Clause 1 fails: no non-canonical dual coexists.  So
the spatial moment rides as a property — a ``spatial_moments`` axis composed
onto BOTH angular field types (``_compose_spatial_moments`` in the bulk; the
flat face-buffer moment tail minted by
:attr:`~orpheus.sn.mesh.augmented_mesh.SNMesh.boundary_face_layout` on the boundary) —
rather than as its own field type.  A ``BoundaryMomentField`` leaf whose
partner-check added nothing beyond class identity would be the vacuous naming
leaf the criterion warns against; the transverse boundary moment is therefore a
PROPERTY of the boundary field, the call S9 made.

The defer-with-trigger decision (#263)
--------------------------------------

The first-class spatial ``SpatialMomentField`` type is DEFERRED, with an
explicit trigger, under Issue #263.  The trigger is the arrival of a
**non-canonical spatial dual**: a nodal / point-value within-cell (or
face-current) representation enters production AND a modelled, applied
nodal↔modal morphism is written between it and the existing
:class:`~orpheus.numerics.spaces.spatial_moment_space.SpatialMomentSpace`.  Two
concrete arrivals would supply it:

* **Nodal discontinuous Galerkin SN** — nodal Lagrange point-values coexist
  with the modal Legendre coefficients, bridged by the applied Vandermonde
  matrix (truncation content, adjoint, in the algebra).  The
  Hesthaven–Warburton nodal-DG construction is the canonical instance.
* **Nodal diffusion (NEM / SANM / ANM)** — transverse-Legendre moments coexist
  with face partial currents, bridged by the coupling coefficients (a modelled
  morphism).  This is the strongest spatial for-case, though it is not on the
  current roadmap.

Until such a dual exists, every order-expansion in the codebase is a
PARAMETER within one tower (the ``per_axis = 2 → 3`` widening for higher-order
SN, hierarchical-Legendre / hp-FEM degree :math:`p`, p-multigrid — all single
hierarchical towers where prolong is inclusion and restrict is its adjoint,
WITHIN one representation), not a new type.  When the dual arrives, the right
move is to lift ``SpatialMomentField`` and its nodal dual into the
:class:`~orpheus.transport.fields._bases.MomentField` family ABC, mirroring the
:math:`M` / :math:`R` pair — the ABC already exists as a thin family marker
anticipating exactly this second instance.  p-adaptivity that needs modelled
``prolong`` / ``restrict`` operators flips toward typed OPERATORS within one
family (like ``truncate``), not a new field type, because those morphisms are
canonical within one Legendre tower.


.. _bc-tensor-primitives:

Boundary conditions as Wave-0 / Wave-1 primitives
=================================================

The boundary trace-law refactor (12-wave effort,
``.claude/plans/transient-giggling-cake.md``, archival narrative
in :ref:`theory-boundary-conditions`) lifted boundary-condition
**realisation** into the same operator algebra documented in this
page: every *realised* BC IS a Wave-0
:class:`~orpheus.numerics.operator.LinearOperator`, composable
with the streaming / scattering / fission operators through the
same algebra dunders.

.. note::

   **A boundary law is NOT itself a Wave-0 operator.** Post Issue
   #186 / B3 + β2 (2026-05-11), the
   :class:`~orpheus.geometry.boundary.BoundaryTraceLaw` ABC is a
   **pure descriptor** carrying no :meth:`apply` method (the
   :class:`LinearOperator` inheritance was dropped). The
   table below lists the **realised** output produced by
   :class:`~orpheus.sn.boundary.realizer.SNBoundaryRealizer` — that
   output IS a Wave-0 :class:`LinearOperator`. The law-to-operator
   transition is the §16A.3 three-layer architecture's realiser
   bridge, enforced statically by the type system rather than by
   convention. See :ref:`bc-trace-law-descriptor-model` for the
   design rationale.

Three new primitives in :mod:`orpheus.numerics.operator`
(Wave 0 of the BC refactor) plus one SN-specific primitive in
:mod:`orpheus.sn.boundary.angular` (Wave 1) form the realisation
target set:

.. list-table:: BC realization primitives — the §15.2 G_α geometric operators
   :header-rows: 1
   :widths: 26 36 38

   * - Primitive
     - Mathematical action
     - Realizes which BC
   * - :class:`~orpheus.numerics.operator.PermutationOperator`
     - ``np.take(x, perm, axis=axis)`` with optional
       involution-detection
     - :class:`~orpheus.geometry.boundary.ReflectiveBoundary`
       (specular reflection via
       ``quadrature.reflection_index(axis)``)
   * - :class:`~orpheus.numerics.operator.IncomingOrdinateMaskTensor`
     - Sparse mask zeroing only the inflow ordinates at one face;
       preserves outflow rows. Self-adjoint + idempotent
       (projector).
     - :class:`~orpheus.geometry.boundary.VacuumInflow` (§16A.10
       trace-correct vacuum)
   * - :class:`~orpheus.numerics.operator.PeriodicWrapOperator`
     - Today: angular identity (the SN realization of
       :class:`~orpheus.geometry.boundary.PeriodicBoundary`).
       Reserved for future spatial-pushforward extension.
     - :class:`~orpheus.geometry.boundary.PeriodicBoundary`
   * - :class:`~orpheus.sn.boundary.angular.AngularAverageOperator`
     - Cosine-weighted Lambertian average over an outgoing
       hemisphere: scalar-broadcasts to all inflow ordinates.
       SN-specific (depends on
       :class:`~orpheus.numerics.quadrature.Quadrature`).
     - :class:`~orpheus.geometry.boundary.WhiteBoundary`
   * - :class:`~orpheus.sn.boundary.angular.IncomingSourceOperator`
     - Ignores the outgoing flux; returns the prescribed source
       value via :meth:`InflowSourceSpec.evaluate`.
     - :class:`~orpheus.geometry.boundary.PrescribedInflow`

The rank-N case (Marshak / partial-current) is built directly
through the algebra dunders:

.. math::
   :label: bc-rank-n-as-sum-of-products

   R \;=\; \sum_\alpha G_\alpha \otimes A_\alpha,
   \qquad
   G_\alpha \in \{\text{permutation, mask, average, wrap, identity}\},
   \quad A_\alpha \in [0, 1].

.. vv-status: bc-rank-n-as-sum-of-products documented

For the SN realisation, the tensor product :math:`G_\alpha \otimes
A_\alpha` collapses to a scalar amplification:
:class:`~orpheus.numerics.operator.ScaledOperator` ``(α, G)``. The
canonical Marshak BC

.. math::

   R_{\text{Marshak}}
   \;=\; c_1 \, G_{\text{refl}} \;+\; c_2 \, G_{\text{diff}}
   \qquad (c_1 + c_2 \leq 1)

becomes the Wave-0 expression
``c1 * PermutationOperator(perm) + c2 * AngularAverageOperator(...)``
— an :class:`~orpheus.numerics.operator.OperatorSum` of
:class:`~orpheus.numerics.operator.ScaledOperator`-wrapped
primitives — **after** the descriptor tree has been realised.

.. _bc-descriptor-tree-vs-operator-tree:

The descriptor-tree algebra is a separate type family
-----------------------------------------------------

The rank-N composition is **not** built directly on the operator
tree. Two distinct algebras are layered:

.. list-table:: Two algebras, two type families
   :header-rows: 1
   :widths: 22 30 24 24

   * - Layer
     - Node types
     - Has ``apply``?
     - When it's built
   * - Descriptor tree
     - :class:`~orpheus.geometry.boundary.BoundaryTraceLaw` leaves;
       :class:`~orpheus.geometry.boundary.LawSum` /
       :class:`~orpheus.geometry.boundary.LawScaled` composers
     - **No** (static type error to call)
     - User code at law-declaration time
   * - Operator tree
     - :class:`~orpheus.numerics.operator.LinearOperator` leaves;
       :class:`~orpheus.numerics.operator.OperatorSum` /
       :class:`~orpheus.numerics.operator.ScaledOperator` composers
     - **Yes** (1-arg :meth:`apply`)
     - Realiser code at face-resolution time

The user writes ``0.3 * spec + 0.7 * white`` — a descriptor tree.
:func:`~orpheus.geometry.boundary.realize_recursively` (method-blind
since #290 P7b) walks the tree, realises each leaf via the
explicitly-passed realizer (e.g.
:class:`~orpheus.sn.boundary.realizer.SNBoundaryRealizer`), and
re-assembles the result through the matching Wave-0 composers
(:class:`OperatorSum` ↔ :class:`LawSum`,
:class:`ScaledOperator` ↔ :class:`LawScaled`). The output is the
operator tree.

**Mixing the two algebras is a type error.** ``LawSum + OperatorSum``
or ``LinearOperator + LawScaled`` are rejected statically (no
matching dunder) — callers MUST realise the descriptor tree first
before composing with the operator algebra. This separation is
the load-bearing payoff of the Issue #186 / B3 + β2 cleanup:
"this is a law, that is an operator" becomes a type-system claim
rather than a convention.

The pre-refactor implementations had two predecessors that
converged on this design:

* **Wave 11** retired the dedicated
  ``MixedBoundaryOperator(components: list[tuple[float, BC]])``
  class — its delayed-realisation pattern broke down once vacuum
  needed per-face inflow indices that the bare-law container
  could not deliver.
* **β1 interim** (Issue #186 / B3, pre-cleanup) kept
  :class:`LinearOperator` on :class:`BoundaryTraceLaw`, so
  ``0.3 * spec + 0.7 * white`` produced an :class:`OperatorSum`
  with raw-law leaves. β1 was algebraically equivalent to β2 but
  conflated the two type families — the type checker could not
  distinguish a not-yet-realised "operator" from a real operator.
  β2 separates them explicitly (the present design).

See :ref:`bc-realize-recursively` for the walker semantics and
:ref:`bc-rank-n-algebra` for the rank-N algebra in detail.


.. _wave-t-tensor-network:

Tensor-Network Decomposition of SN Operators (Wave T)
=====================================================

Wave T (May 2026, commits ``fa13e78`` / ``0b2848b`` / ``9f85c5d`` /
``03bcdba`` / ``cb18fdb`` / ``c55b505`` / ``90e7d4e``) lifted the four
SN operator leaves — boundary realizers, fission, scattering,
streaming — from procedural single-axis numpy bodies into the
operator-algebra types documented above
(:class:`~orpheus.numerics.operator.TensorProductOperator`,
:class:`~orpheus.numerics.operator.SumOfTensorProductsOperator`,
:class:`~orpheus.numerics.operator.OperatorSum`). The migration is
the consumer side of the Wave-0 + Depth B D-B infrastructure
(:class:`~orpheus.numerics.operator.TensorProductOperator` shipped at
commit ``bc1253e``, 2026-05-10;
:class:`~orpheus.numerics.space.TensorProductSpace` shipped at commit
``c2f968a``, 2026-05-27; the first production consumer was the D-B+1
specular BC at ``boundary_realizer.py:164-166``).

What landed is **not** the uniform :math:`A = \sum_k A_x^{(k)} \otimes
A_\omega^{(k)} \otimes A_g^{(k)}` aspiration of Grand Report v3
§15-§16A. The shipped state is **five algebraically distinct shapes**,
chosen per operator by what the underlying physics actually couples.
Future agents who assume "all SN operators are
:class:`SumOfTensorProductsOperator`" — including Wave O
(`Issue #208 <https://github.com/deOliveira-R/ORPHEUS/issues/208>`_)
operator-role typing — will be wrong. This section names the master
condition that decides the shape, catalogues which shape each operator
uses, and documents the architectural rationale for the per-direction
streaming split.


Key Facts (Wave T)
------------------

- **The SN flux state lives on a tensor-product space** :math:`V = X
  \otimes \Omega \otimes G` (Grand Report v3 §15 line 2003-2019). The
  shipped array layout ``(N, ng, nx, ny)`` is the implicit numpy
  realisation: the angular axis :math:`\Omega` is leading, the group
  axis :math:`G` is next, and the spatial axis :math:`X` trails (see
  :ref:`theory-sn-index-convention`).

- **Five algebraic shapes** are now in production simultaneously,
  selected by the per-operator physics coupling — see
  :ref:`wave-t-shape-table` below for the catalogue.

- **The MA-Q1 master condition** (load-bearing for every future
  consumer):

  .. epigraph::

     :class:`SumOfTensorProductsOperator` (SOTP) requires Cartesian-
     product per-axis decomposition: every summand factors as a
     product of independent per-axis operations. *Coupled physics* —
     per-material XS lookup that ties group to spatial cell,
     sequential WDD recurrence that ties spatial cells, M-M half-grid
     recurrence that ties angular ordinates — falls back to
     :class:`OperatorSum` over bespoke :class:`LinearOperator`
     summands, **NOT** SOTP.

- **Zero production consumers** of
  :class:`SumOfTensorProductsOperator`. The §15.2 SOTP form is
  contradicted by the actual coupling structure of scattering (T.3)
  and streaming (T.4). Only T.1 (BC realizers) and T.2 (fission rank-1)
  cleanly admit the clean tensor-product factorisation.

- **Wave O typing constraint**: operator-role types
  (``BulkOperator`` / ``FullOperator`` / ``BoundaryOperator``
  Protocols, Issue #208) **MUST** accept non-SOTP summands. Any
  contract that requires "all summands are
  :class:`TensorProductOperator`" forecloses scattering and streaming.

- **In-sweep streaming (the retired per-direction split, #238)**: the
  forward (μ_x > 0) and backward (μ_x < 0) WDD recurrences are walked in
  ONE bidirectional pass inside
  :meth:`~orpheus.sn.loss_representation._OneDimScanWalk._apply_walk` (the
  fused ``(L+C)ψ`` matvec, the apply-direction twin of the sweep). Wave T
  briefly exposed the two directions as separately-applicable typed leaves
  (``M_spatial = _SpatialSweepDirection(+1) + _SpatialSweepDirection(-1)``)
  to anticipate Wave-O adjoint / DSA-class consumers — but #240's adjoint
  uses the fused ``loss_action_transpose`` and no production code ever
  applied the leaves separately, so #238 retired the split. The
  per-direction coupling structure it documented (the backward sweep's seed
  depends on the forward sweep's outer-face WDD outflow) is still real and
  is why the matvec is one shared pass, not two independent forward sweeps.

- **Hybrid 2-D Cartesian**: T.4 lifted 1-D only. The 2-D Cartesian
  path (then ``StreamingOperator._apply_2d_cartesian``) stayed
  procedural FD with cell-centre-proxy semantics, guarded by a
  defensive source-hash pin (A2D-1) against silent author drift.
  *(Superseded: O.4b replaced the FD path with the trace-correct
  graph walk; S6.3 moved the walk onto the loss representation; at
  S6.4(a) the A2D-1 pin was RETIRED — its tripwire job transferred to
  the ``window ≡ full`` matvec output oracle in
  ``test_2d_full_field_oracle.py``, which catches actual drift by
  ``assert_array_equal`` against the structurally-distinct full-field
  oracle instead of pinning source text on what is now a SHARED
  octant walk.)*

.. vv-status: wave-t-shape-table documented


.. _wave-t-shape-table:

Per-operator shape catalogue
----------------------------

The five algebraic shapes that ship today, grouped by Wave T substep.
Each row names the operator, the algebraic shape its kernel/apply
takes, a concrete example, and the physics coupling that forced the
shape choice.

.. list-table:: Wave T per-operator shape catalogue
   :header-rows: 1
   :widths: 18 22 30 30

   * - Operator
     - Algebraic shape
     - Example
     - Why this shape
   * - BC realizers (vacuum,
       specular, white,
       albedo, periodic)
     - :class:`TensorProductOperator`
       (single TP)
     - ``IncomingOrdinateMaskTensor(...) &
       IdentityOperator()`` for vacuum;
       ``PermutationOperator(perm, axis=0) &
       IdentityOperator()`` for specular
     - Each BC acts on the ordinate axis; the
       trailing group / face axes broadcast. Per
       §16A.10 ``B = G_patch ⊗ K_omega ⊗ K_g``
       with two factors degenerate to
       :class:`IdentityOperator`.
   * - Fission (:math:`F`)
     - :class:`TensorProductOperator`
       (single rank-1 dyad)
     - ``outer(χ, ReactionRateFunctional(νΣ_f)) &
       IdentityOperator()``
       (:attr:`FissionOperator.kernel`)
     - Per Grand Report v3 §15 line 2008
       :math:`F = |\chi\rangle\langle\nu\Sigma_f|`. The
       group-axis contraction-then-broadcast is
       exactly :class:`RankOneOperator`; spatial
       axes broadcast.
   * - Scattering kernel
       (:math:`S_{\rm aniso}`)
     - :class:`OperatorSum` of
       per-ℓ bespoke leaves
     - ``reduce(add, kernel_summands)`` where each
       summand is
       :class:`_PerLegendreOrderScattering(ell=ℓ)`
       (:attr:`ScatteringOperator.kernel`)
     - **MA-Q1 fallback**: the per-material per-ℓ
       einsum
       :meth:`MaterialXSField.apply_legendre_scattering_moments`
       couples the group axis (matrix multiply on
       :math:`\Sigma_{s,\ell}[g'\to g]`) with the
       spatial axis (via
       :attr:`cells_by_material` indexing). No
       SOTP factorisation respects disjoint axes;
       the §15.2 SOTP target form fails the
       :class:`TensorProductOperator` contract.
   * - Streaming spatial
       part (in-sweep;
       #238 retired the
       ``M_spatial`` leaf)
     - Fused WDD recurrence in
       :meth:`~orpheus.sn.loss_representation._OneDimScanWalk._apply_walk`
     - The forward (μ_x > 0) + backward (μ_x < 0)
       sweeps walked in ONE bidirectional pass
       (formerly the ``M_spatial`` per-direction
       :class:`OperatorSum`)
     - **MA-Q1 fallback**: the WDD recurrence
       :math:`\psi_{\text{face,out}} = 2\bar\psi
       - \psi_{\text{face,in}}` sequentially
       couples cells along x. It is NOT a clean
       :math:`(D_x \otimes \Omega_x \otimes I_g)`
       3-factor TP — the sweep operator is the
       leaf factor. The two directions cannot be
       applied independently (the backward seed
       depends on the forward outer-face outflow),
       which is why the matvec is one shared pass.
   * - Streaming angular
       redistribution
       (in-sweep; #238
       retired the
       ``M_angular_redist``
       leaf)
     - In-sweep Morel–Montry thread inside
       :meth:`~orpheus.sn.loss_representation._OneDimScanWalk._apply_walk`
       (sphere / cylinder; zero for slab / 2-D
       Cartesian)
     - The per-cell M-M contribution from
       :meth:`PoleAngularClosureBase.cell_contribution`,
       added to the cell balance during the walk
     - **MA-Q1 fallback**: the M-M half-grid
       recurrence (Hébert 2009 §3.9.4
       Eqs. 3.432-3.435) sequentially couples
       angular ordinates ``α_{m+1/2}`` from
       ``α_{m-1/2}`` with σ_t-dependent
       absorption coefficients. Not a diagonal
       angular factor; a 3-factor TP wrap would
       false-assert separability the recurrence
       doesn't support. Verified end-to-end by the
       anisotropic curvilinear MMS
       (:ref:`sn-mms-curvilinear-aniso-verification`).


The MA-Q1 master condition
--------------------------

The pattern across T.3, T.4-spatial, and T.4-curvilinear is the same:
*coupled physics produces summands that fail the disjoint-axes
contract of* :class:`TensorProductOperator`. Naming this explicitly
prevents future agents from re-attempting the §15.2 SOTP form on each
of these operators.

.. (vv-status rationale) Master condition gate for Wave-O typing
   decisions. Verified by the absence of SOTP-shaped consumers in
   production after T.3 + T.4 land — exhaustively documented in
   ``.claude/plans/wave_t_tensor_network.md`` §6 T.3 + T.4 deviations.
.. vv-status: wave-t-ma-q1-master-condition documented

.. math::
   :label: wave-t-ma-q1-master-condition

   \text{SOTP applies} \;\Longleftrightarrow\;
   \text{each summand factors as} \;
   f(x_1,\dots,x_d) \;=\; f_1(x_1)\otimes\cdots\otimes f_d(x_d).

When the physics violates the right-hand side — and three of the four
Wave-T-touched operators do — the algebraic home is
:class:`OperatorSum` over :class:`LinearOperator` summands, NOT
:class:`SumOfTensorProductsOperator`. The §15.2 target form is
*aspirational* in the grand report; Wave T documents that the actual
coupling structure of multigroup transport with per-material
cross-sections does not admit it for scattering and streaming.

**Three coupled-physics archetypes** ship in Wave T:

1. **Per-material XS coupling** (T.3 scattering). The per-material
   einsum :math:`\sum_{g'}\Sigma_{s,\ell}^{m(\vec r)}[g'\to g]
   \phi_{\ell,g'}^{m}(\vec r)` ties the group axis (matrix multiply)
   to the spatial axis (per-cell material id lookup). The factor
   :math:`\Lambda_\ell` cannot be written as a group-axis-only
   operator without information loss.

2. **Sequential WDD recurrence** (T.4 streaming, spatial). The
   diamond-difference closure :math:`\psi_{\text{face,out}} =
   2\bar\psi_{\text{cell}} - \psi_{\text{face,in}}` makes the cell
   :math:`i+1` value depend on the cell :math:`i` value. A
   per-direction sweep summand IS the WDD recurrence as a single
   :class:`LinearOperator` — not a factor on a per-cell tensor axis.

3. **M-M half-grid recurrence** (T.4 streaming, curvilinear angular).
   The Carlson-Morel-Montry α-coefficients (Hébert 2009 §3.9.4
   Eqs. 3.432-3.435) recur sequentially along the angular axis within
   each μ-level: :math:`\alpha_{m+1/2}` depends on
   :math:`\alpha_{m-1/2}` and on σ_t. The leaf factor is the entire
   recurrence — a single :class:`LinearOperator` — not a diagonal
   angular operator.

In each case, the algebraic home is the SAME — :class:`OperatorSum`
over bespoke :class:`LinearOperator` summands — and the
:class:`TensorProductOperator` form is structurally inaccessible
without information loss.


.. _wave-t-streaming-deep-dive:

Streaming deep dive — in-sweep WDD recurrence (the retired per-direction split)
-------------------------------------------------------------------------------

.. note:: #238 — the ``M_spatial`` / ``M_angular_redist`` typed-leaf split was retired.

   Wave T briefly exposed the streaming matvec as separately-applicable
   typed leaves: ``StreamingOperator.M_spatial`` (an :class:`OperatorSum`
   of two per-direction-sign sweep summands) and
   ``StreamingOperator.M_angular_redist`` (a bespoke curvilinear
   angular-redistribution leaf). The split was designed to anticipate
   Wave-O adjoint propagation, DSA-class preconditioners, and
   per-direction debugging. **#238 retired it**: no production code ever
   applied the leaves separately (the #240 G-adjoint rides the fused
   ``loss_action_transpose``; the open #200 block-inverse preconditioner
   and #2 DSA never landed), so keeping the leaves alive solely to feed
   their own structural tests was the same orphan smell one level down.
   The streaming + curvilinear Morel–Montry angular redistribution is now
   computed **in-sweep** inside the fused matvec
   :meth:`~orpheus.sn.loss_representation._OneDimScanWalk._apply_walk` (the
   apply-direction twin of the sweep), and the angular-redistribution term
   is verified end-to-end by the anisotropic curvilinear MMS
   (:ref:`sn-mms-curvilinear-aniso-verification`,
   ``catches("ERR-026")``) — the surviving structural-independence ground.

The *physics* the split documented is unchanged and load-bearing: the
WDD spatial recurrence and the M-M angular recurrence are both
sequentially coupled, which is why neither admits a clean
:class:`TensorProductOperator` factorisation (the MA-Q1 master condition
above). The two equations below pin that coupling, and they are the
reason the matvec is a single sequential walk rather than a stack of
independent per-axis broadcasts.

**Where the WDD recurrence comes from**. For a single ordinate with
:math:`\mu_x > 0`, the discrete cell balance over cell :math:`i`
(width :math:`\Delta x_i`, total cross-section :math:`\Sigma_t(i)`,
cell-averaged in-group source :math:`\bar Q_i`) is

.. math::

   |\mu|\,\bigl[(\psi_{\text{face,out}})_i - (\psi_{\text{face,in}})_i\bigr]
   \;+\; \Sigma_t(i)\,\Delta x_i\,\bar\psi_i
   \;=\; \Delta x_i\,\bar Q_i .

This single equation carries **two** unknowns —
:math:`(\psi_{\text{face,out}})_i` and the cell average
:math:`\bar\psi_i`. The diamond-difference (DD) closure supplies the
second relation, asserting the cell average is the arithmetic mean of
the two faces,

.. math::

   \bar\psi_i \;=\;
     \tfrac12\bigl[(\psi_{\text{face,in}})_i + (\psi_{\text{face,out}})_i\bigr]
   \;\Longleftrightarrow\;
   (\psi_{\text{face,out}})_i \;=\; 2\,\bar\psi_i - (\psi_{\text{face,in}})_i .

Substituting the closure into the balance and solving for the cell
average gives the **forward WDD recurrence**:

.. math::
   :label: wdd-forward-recurrence

   \bar\psi_i \;=\;
     \frac{\Delta x_i\,\bar Q_i + |\mu|\,(\psi_{\text{face,in}})_i}
          {\Delta x_i\,\Sigma_t(i) + |\mu|},
   \qquad
   (\psi_{\text{face,out}})_i \;=\;
     2\,\bar\psi_i \;-\; (\psi_{\text{face,in}})_i

with :math:`(\psi_{\text{face,in}})_{i+1} =
(\psi_{\text{face,out}})_i` (the outflow of cell :math:`i` is the
inflow of cell :math:`i+1`).

**Why this forbids a tensor product**. The recurrence is a
*sequential* dependence: the cell-:math:`i+1` average cannot be
formed until the cell-:math:`i` outflow is known, which itself
depends on cell :math:`i-1`, and so on back to the boundary inflow.
Unrolling the recurrence to its closed form makes the structure
explicit — the cell average is a lower-triangular linear functional
of the upstream source and the inflow face:

.. math::

   \bar\psi_i \;=\;
     \frac{|\mu|}{\Delta x_i\Sigma_t(i)+|\mu|}\,\psi_{\text{bdy,in}}
     \prod_{j<i}\frac{2|\mu| - \Delta x_j\Sigma_t(j) - |\mu|}
                     {\Delta x_j\Sigma_t(j)+|\mu|}
     \;+\;\bigl(\text{source terms}\bigr),

i.e. the action on the spatial axis is the *whole* lower-triangular
sweep operator :math:`T_x`, not a per-cell diagonal that broadcasts.
There is no factorisation :math:`(D_x \otimes \Omega_x \otimes I_g)`
in which :math:`D_x` acts independently on each spatial cell: the
off-diagonal products :math:`\prod_{j<i}(\cdots)` carry information
*between* cells and *depend on the ordinate* :math:`\mu` through the
denominators. The spatial factor and the angular index are entangled
inside the recurrence — the disjoint-axes contract of
:class:`TensorProductOperator` (:eq:`wave-t-ma-q1-master-condition`)
fails. The cell-balance algebra at
:func:`orpheus.transport.spatial.cell_balance.cell_balance_for_streaming`
hides this recurrence inside the named denom-numer primitives, but
the recurrence is the load-bearing structure that makes the streaming
matvec a sweep rather than a broadcast.

.. vv-status: wdd-forward-recurrence documented

**Forward-backward coupling at the outer face**. The forward (μ_x > 0)
and backward (μ_x < 0) sweeps cannot be applied independently: the
backward sweep's seed depends on the forward sweep's outer-face WDD
outflow. Concretely, the forward sweep marches from the inner boundary
to the outer boundary, terminating with the outer-face outflow
:math:`(\psi_{\text{face,out}})_{N-1}`; the backward sweep then marches
from the outer boundary back inward, and on a reflective or curvilinear
mesh its seed inflow at the outer face is determined by that same
forward outflow. The two directions therefore cannot run as two parallel
independent forward sweeps — the data dependency from the forward
terminus into the backward seed is intrinsic to the recurrence. This is
why
:meth:`~orpheus.sn.loss_representation._OneDimScanWalk._apply_walk`
walks both directions in ONE bidirectional pass, regardless of whether
the directions are exposed as separate operator leaves. (Exposing them
as leaves does not remove the coupling — it only forces the leaves to
share state, which is precisely the smell that motivated the
now-retired orchestrator; see :ref:`wave-t-orchestrated-apply`.)

.. note::

   The forward-outflow-seeds-backward-inflow coupling was the
   *pre*-O.4a.2 mechanism, where the backward sweep's inflow was
   ``bc_outer.apply(forward_outflow)`` inside one matvec. Wave O step
   O.4a.2 (Issue #208) **deleted that intra-call reflective re-apply** —
   the boundary law :math:`B` is now a sibling :math:`-B` operator and
   the backward sweep reads the *given* outer inflow trace directly. The
   forward outflow today feeds the **outflow self-consistency defect** on
   the outflow trace row, not a reflected inflow seed.
   See :ref:`bc-extraction`.


.. _wave-t-orchestrated-apply:

One bidirectional pass — design rationale and the retired per-direction split
-----------------------------------------------------------------------------

The fused matvec
:meth:`~orpheus.sn.loss_representation._OneDimScanWalk._apply_walk` runs
the bidirectional sweep **once**, returning the full :math:`(L+C)\,\psi`.
This is the whole production story today; the rest of this section records
*why* Wave T briefly exposed a richer surface, *what* that surface looked
like, and *why* #238 retired it. The history is load-bearing: it
documents a structural fact — that the per-direction split could never
have been cheaper or simpler than the fused pass — so a future session
does not re-introduce the split on the mistaken belief that it buys
modularity for free.

What Wave T originally built (Design B)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Wave T exposed the streaming matvec as a **subclass of**
:class:`OperatorSum`, ``_MSpatialOperatorSum``, whose two summands were
``_SpatialSweepDirection(+1)`` (the μ_x > 0 forward sweep) and
``_SpatialSweepDirection(-1)`` (the μ_x < 0 backward sweep), so that

.. math::

   M_{\rm spatial} \;=\;
     \texttt{_SpatialSweepDirection}(+1)
     \;+\;
     \texttt{_SpatialSweepDirection}(-1) .

Carrying the two directions as named ``.a`` / ``.b`` attributes of an
:class:`OperatorSum` made them visible to type introspection by Wave O
(:ref:`bc-extraction`), adjoint propagation, and any DSA-class
preconditioner that might want to address one direction at a time. The
subclass **overrode** :meth:`OperatorSum.apply`, because the default
implementation

.. math::

   \texttt{OperatorSum.apply}(x) \;=\;
     \texttt{self.a.apply}(x) \;+\; \texttt{self.b.apply}(x)

would have cost **1.5× the unified matvec walltime**. The reason is the
forward-backward coupling of :ref:`wave-t-streaming-deep-dive`: each
standalone ``_SpatialSweepDirection.apply`` internally ran the *entire*
bidirectional sweep (it had to, to obtain the seed coupling) and then
masked the opposite-direction ordinates to zero. Summing the two
standalone summands therefore re-ran the forward sweep twice. The
orchestrator's override ran the bidirectional sweep once via the shared
``_MSpatialOperatorSum._compute_LpC`` and returned the full
:math:`(L+C)\,\psi`, avoiding the duplication; the standalone
per-direction summands were preserved only as a slow fallback for
testing, adjoint inspection, and per-direction debugging.

The forward-sweep outer-face WDD outflow that had been a hidden local
of the legacy unified matvec was lifted by the orchestrator into the
*named shared state* of ``_compute_LpC`` (``coding-elegance`` Pattern 6
— single source of truth: hidden coupling points must become named).
The full bidirectional matvec is mathematically equivalent to
:math:`M_{x,+}\,\psi + M_{x,-}\,\psi`; the orchestrator returned that
value bit-exact (preserving the unified matvec's reduction order), while
the masked per-direction summands summed to the same value at
FP-non-associativity ULP.

The five anticipated consumers — and why none materialised
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The split was speculative architecture: it was built to *anticipate*
future consumers that would want a separately-applicable per-direction
or per-term operator leaf. The honest record of what was expected, and
what actually happened, is:

.. list-table:: Anticipated consumers of the per-direction / per-term split
   :header-rows: 1
   :widths: 26 38 36

   * - Anticipated consumer
     - Why a separate leaf was expected to help
     - What actually shipped (why it didn't materialise)
   * - **Wave-O adjoint** (Issue #208 / #240 G-adjoint)
     - A per-direction leaf was expected to make
       :math:`M_{\rm spatial}^{\mathsf T}` introspectable
       direction-by-direction, so the adjoint could be assembled by
       transposing each summand.
     - The #240 G-adjoint rides the **fused**
       :meth:`~orpheus.sn.loss_representation._OneDimScanWalk.loss_action_transpose`
       — the apply-direction twin transposed as a whole walk. The
       transpose of a lower-triangular sweep is an upper-triangular
       sweep over the *same* coupling; splitting it per direction first
       gains nothing and re-introduces the shared-state coupling.
   * - **#2 consistent DSA**
     - A diffusion-synthetic accelerator was expected to consume an
       isolated spatial-streaming leaf as the operator to precondition.
     - DSA consumes the **fused residual** of :math:`(L+C-S)` plus a
       *separate* diffusion operator :math:`A_{\rm diff}` built
       in-algebra; it never needs the streaming term split from
       collision, let alone the forward direction split from the
       backward. That diffusion operator now **exists** —
       :math:`A_{\rm diff} = L + C - S - B` on the scalar composite
       (:mod:`orpheus.diffusion.operators`, #290 P4;
       :doc:`/theory/diffusion_1d`), inverted by an explicit
       :class:`~orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator`
       and correctness-safe by construction (its low-order correction
       :math:`\to 0` at convergence). The accelerator (Issue #2) itself
       is still unbuilt, but the architecture it decided on is this
       in-algebra diffusion operator, not a per-direction streaming
       leaf.
   * - **#200 block-inverse preconditioner**
     - A per-direction leaf was expected to feed a block-inverse Krylov
       preconditioner that addressed each direction block.
     - The full Morel–Montry sweep is *already* the natural
       :math:`O(N)` exact inverse :math:`(L+C)^{-1}` (the sweep solves
       in one pass). A spatial/angular block split of an operator whose
       inverse is already a single cheap pass would be **weaker, not
       cheaper** — a block preconditioner approximates an inverse that
       the sweep computes exactly. #200 remains open and, when it lands,
       has no reason to want the split.
   * - **Per-direction debugging**
     - Inspecting one sweep direction in isolation while debugging.
     - The fused walk is debuggable directly (the per-level / per-cell
       visits are observable inside the single pass); a standalone
       direction summand that re-runs the whole bidirectional sweep and
       masks the other direction is a *worse* debugging surface than the
       single pass, because the masking hides the coupling under test.
   * - **Slow per-direction test fallback**
     - A standalone per-direction ``apply`` as a structural cross-check
       on the fused pass.
     - The cross-check that mattered — that the fused
       :math:`(L+C)\,\psi` is correct — is supplied by the
       anisotropic curvilinear MMS
       (:ref:`sn-mms-curvilinear-aniso-verification`), a
       *structurally-independent* L1 ground. A bit-identity invariant
       between the fused pass and the sum of its own per-direction
       summands (see :ref:`wave-t-curvilinear-deep-dive`) only verified
       that the split *reconstructed itself*, not that either branch was
       correct.

Why #238 retired the split — the orphan-smell one level down
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After #206 moved the 1-D matvec walk into the loss representation, the
``_SpatialSweepDirection`` / ``_MSpatialOperatorSum`` /
``AngularRedistributionOperator`` leaves had **zero production
consumers**. Every production matvec went through the fused walk; the
leaves existed only to be applied by their own structural tests — the
algebra-decomposition invariant
:math:`(L+C)\,\psi \equiv M_{\rm spatial}\,\psi + M_{\rm angular\_redist}\,\psi`
(see :ref:`wave-t-curvilinear-deep-dive`). This is the classic *orphan
smell*: machinery kept alive solely to feed the tests that exist solely
to verify that machinery. Cardinal Rule 2 (architecture is critical)
treats a self-referential test loop as a code smell — the tests prove
the decomposition is internally consistent without proving anything a
production consumer relies on.

The decisive observation is that the split was never going to be the
cheaper *or* the more modular path, for a structural reason rather than
an incidental one:

- **It could not be cheaper.** The forward-backward coupling
  (:ref:`wave-t-streaming-deep-dive`) means a standalone per-direction
  apply must run the whole bidirectional sweep anyway. The orchestrator
  existed precisely to undo the 1.5× penalty of the naïve sum — i.e. to
  claw back to the cost of the single fused pass it was wrapping. The
  fused pass is the floor; the split can at best tie it.
- **It could not be more modular.** Exposing the directions as leaves did
  not decouple them — it forced them to *share named state*
  (``_compute_LpC``'s lifted outflow), so the "modular" surface was a
  pair of objects that could not be evaluated independently. That is the
  illegal state ``coding-elegance`` Pattern 4 warns against:
  representing two summands as separable when the math makes them
  inseparable.

So #238 removed the orchestrator and both leaves entirely. The single
bidirectional pass survives as the production path; the 1.5× cost the
override avoided is now moot — there is only one pass to run — and the
single-source-of-truth property is satisfied *more directly*: one walk,
one source, no orchestrator coordinating two summands that were never
independent. The angular-redistribution term that the curvilinear leaf
isolated is computed in-sweep (see :ref:`wave-t-curvilinear-deep-dive`)
and verified end-to-end by the anisotropic curvilinear MMS
(:ref:`sn-mms-curvilinear-aniso-verification`, ``catches("ERR-026")``)
— the surviving structural-independence ground.

.. note:: **If a future consumer genuinely needs a per-direction or
   per-term leaf.** Should a #200 block-inverse preconditioner or a #2
   DSA variant ever surface a *real* need for a separately-applicable
   spatial / angular leaf (one that does not re-run the full
   bidirectional sweep), the correct move is **not** to resurrect
   ``_MSpatialOperatorSum``. The structural obstruction above — that the
   forward-backward coupling makes the directions inseparable — has to be
   addressed first: the consumer would need a formulation in which the
   coupling itself is the object being preconditioned (e.g. the
   in-algebra diffusion operator for DSA — now built as
   :math:`A_{\rm diff} = L + C - S - B`, #290 P4,
   :doc:`/theory/diffusion_1d` — per the architecture decided on
   Issue #2), not a re-split of the sweep into directions that secretly
   share state.
   The fused
   :meth:`~orpheus.sn.loss_representation._OneDimScanWalk.loss_action`
   /
   :meth:`~orpheus.sn.loss_representation._OneDimScanWalk.loss_action_transpose`
   pair is the principled surface; build the new consumer against it.


.. _wave-t-curvilinear-deep-dive:

Curvilinear angular redistribution — in-sweep Morel–Montry thread
-----------------------------------------------------------------

For sphere / cylinder geometries the curvilinear M-M (Morel–Montry)
half-grid angular redistribution is woven into the cell balance during
the fused walk (it returns zero for slab / 2-D Cartesian). The per-cell
M-M coefficients come from
:meth:`PoleAngularClosureBase.cell_contribution` (Pattern 6 — single source
of truth for the M-M coefficients). #238 retired the bespoke
``AngularRedistributionOperator`` leaf that re-walked the matvec only to
isolate this term; the redistribution is the same in-sweep computation,
now without a separately-applicable wrapper.

**Where the angular-redistribution term comes from**. In curvilinear
geometry the streaming operator acquires a term with no slab analogue:
the angular derivative that accounts for the rotation of the local
direction frame as a particle streams along a curved trajectory. For
the sphere this is :math:`\frac{1-\mu^2}{r}\,\partial\psi/\partial\mu`;
for the cylinder it is :math:`-\frac1r\,\partial(\xi\psi)/\partial\varphi`.
Discretised over the cell volume (see the canonical derivation in
:doc:`discrete_ordinates`, Step 2–3, Eq. :eq:`balance-general`), it
becomes a *half-grid* difference between angular faces
:math:`m\pm\tfrac12`,

.. math::

   \int_{V_i}\frac{1-\mu^2}{r}\frac{\partial\psi}{\partial\mu}\,dV
   \;\approx\;
   \alpha_{m+\frac12}\,\psi_{m+\frac12}
   \;-\;
   \alpha_{m-\frac12}\,\psi_{m-\frac12},

with the geometry factor :math:`\Delta A_i / w_n` restoring per-ordinate
flat-flux consistency (without it, the cancellation that makes a spatially
uniform angular flux exact holds only in the *sum* over ordinates, not
per ordinate — the Morel–Montry flux dip near :math:`r=0`). The
half-angle fluxes :math:`\psi_{m\pm 1/2}` are not free unknowns: they are
fixed by a closure that ties each half-angle to its neighbour, which is
the source of the sequential coupling below.

**Why not a tensor product**. Per Hébert 2009 §3.9.4, Eqs. 3.432-3.435,
the M-M closure produces an angular recurrence

.. math::
   :label: mm-half-grid-recurrence

   \alpha_{m+1/2} \;=\;
     f(\alpha_{m-1/2},\;\Sigma_t,\;\psi_{m-1/2,\,\text{upstream}})

within each μ-level :math:`p`. The :math:`\alpha_{m\pm 1/2}` are the
Carlson coupled-pole half-angle coefficients, and the recurrence on
angular ordinates is sequential along the half-grid axis with
σ_t-dependent absorption. The factor that produces
:math:`\alpha_{m+1/2}` from :math:`\alpha_{m-1/2}` IS the entire
recurrence; there is no clean per-angular-axis diagonal factor that
respects the disjoint-axes contract. A 3-factor TP wrap would
**false-assert separability** the recurrence doesn't support
(``coding-elegance`` Pattern 4 — do not represent illegal states).

This is the *angular* analogue of the spatial WDD obstruction of
:ref:`wave-t-streaming-deep-dive`: where the WDD recurrence couples
spatial cell :math:`i+1` to cell :math:`i`, the M-M recurrence couples
angular half-face :math:`m+\tfrac12` to :math:`m-\tfrac12`. Both are
lower-triangular sweeps over a single axis with ordinate- and
material-dependent coefficients, and **both run inside the same fused
walk**, sequentially nested: the outer loop is the spatial sweep, and at
each spatial cell the inner M-M thread advances the angular half-grid.
That nesting is precisely what a :class:`SumOfTensorProductsOperator`
cannot express — a sum of per-axis tensor factors is a flat algebraic
form with no notion of one axis's recurrence running *inside* another's.

.. vv-status: mm-half-grid-recurrence documented

**Per-cell algebra**. The walk visits every :math:`(p,\,i)` pair
(μ-level × spatial cell) and calls
:meth:`PoleAngularClosureBase.cell_contribution`. The cell-balance algebra
at :func:`orpheus.transport.spatial.cell_balance.cell_balance_for_streaming`
decomposes additively into three terms:

.. math::
   :label: wave-t-cell-balance-three-terms

   {\rm denom} \;=\;
     {\rm streaming\_denom\_term} \;+\;
     {\rm angular\_denom\_term} \;+\;
     {\rm collision\_denom\_term}

.. math::

   {\rm numer\_upstream} \;=\;
     {\rm spatial\_upstream\_term} \;+\;
     {\rm angular\_numer\_upstream}

The angular-redistribution contribution to the cell balance is

.. math::

   m_{\rm angular\_redist} \;=\;
     \frac{1}{V_i}\bigl[
        {\rm angular\_denom\_term} \cdot \psi_{\rm cell}
      - {\rm angular\_numer\_upstream}
     \bigr]

with :math:`{\rm angular\_denom\_term} = (\Delta A / w)\,c_{\rm out}`
and :math:`{\rm angular\_numer\_upstream} = (\Delta A / w)\,c_{\rm in}\,
\psi_{m-1/2,\,i,\,g}` per the M-M closure (see
:class:`orpheus.sn.spatial.pole_angular_closure.PoleAngularClosureBase`
for the closure data). It is an interior-cell operation that does not
traverse the spatial boundary; only the spatial sweep writes face
residuals.

.. vv-status: wave-t-cell-balance-three-terms documented


Curvilinear bookkeeping (formerly the M_spatial-via-subtraction smell)
----------------------------------------------------------------------

When the retired ``M_spatial`` leaf existed, its curvilinear value was
defined by *subtracting* the angular-redistribution share from the
unified :math:`(L+C)\,\psi`:

.. math::
   :label: wave-t-mspat-curvilinear-subtraction

   M_{\rm spatial}\,\psi \;=\; (L+C)\,\psi \;-\;
                               M_{\rm angular\_redist}\,\psi
   \qquad (\text{curvilinear})

.. vv-status: wave-t-mspat-curvilinear-subtraction documented

#238 retired both leaves, so this subtractive definition no longer ships:
the fused matvec emits :math:`(L+C)\,\psi` directly with the
angular-redistribution term already folded into the cell balance, and the
spatial / angular split is never materialised. The algebra-decomposition
invariant test that bounded the old subtraction
(:math:`(L+C)\,\psi \equiv M_{\rm spatial}\,\psi + M_{\rm angular\_redist}\,\psi`
at principled-equivalence ULP) retired with the leaves; the surviving
guarantee that the angular term is correct is the anisotropic curvilinear
MMS (:ref:`sn-mms-curvilinear-aniso-verification`), which exercises the
full curvilinear :math:`(L+C)` end-to-end.


The 2-D Cartesian hybrid (Q1)
-----------------------------

T.4 lifted the 1-D path only. The 2-D Cartesian path
(:meth:`StreamingOperator._apply_2d_cartesian`) remains procedural
cell-centred upwind FD with **cell-centre-proxy boundary semantics**:
the matvec body reads ``psi.bulk.values[:, :, 0, iy]`` as the outgoing
trace at xmin and the BC's :meth:`apply(outgoing)` fills the
incoming-direction bulk cells. The
:attr:`psi.boundary.face_view` is currently **passive**: its values
do not enter the bulk computation.

The trace-correct face_view formulation (face_view enters the bulk
computation as the boundary trace, with a boundary residual driving
face_view ↔ bulk consistency) caused a 10% k_inf drift in
experiments — see the existing comment at
``orpheus/sn/operator.py:843-868``. That rewire requires the BC
realizers to gain a "proper composable algebra" — a payload distinct
from T.4's per-direction lift. Bundling the two would violate the
``feedback_unify_after_two_instances`` discipline (only one working
2-D path exists; unify after ≥2 working instances).

**Defensive A2D-1 source-hash pin**. T.4d added a structural
regression test that recorded the source-code signature of
``_apply_2d_cartesian`` and asserted it remained unchanged, so an
accidental modernisation of the 2-D path could not ship silently.
*(RETIRED at S6.4(a): the body became the SHARED ``_OctantWalk``
apply frame, where a source-hash trips on every legitimate refactor
with no behavior signal; the ``window ≡ full`` matvec
``assert_array_equal`` oracle is the successor tripwire.)*


Verification approach
---------------------

Wave T's verification chain combines three independent grounds:

1. **Pre-T.4 snapshot bit-identity** (Route A). The substep T.4a
   captured the value of :meth:`apply` and :meth:`solve` on fixed
   :math:`(\text{seed}, \text{mesh}, \text{material})` triples across
   slab, sphere, cylinder, and 2-D Cartesian, plus 1G / 2G / asymmetric
   :math:`\Sigma_s` / vacuum / white / specular variants. Each
   subsequent T.4 substep is gated on :func:`np.array_equal` against
   those snapshots — the existing numerics are the local
   bit-identity reference.

2. **Principled-equivalence ULP** for cases where reductions reorder.
   Per the ``vv-principles`` skill §"Bit-identity
   vs principled-equivalence": when the operator-algebra fold inserts
   a :func:`np.add` at a different position than the legacy fused
   einsum, the new value is verified by the three-criteria gate
   (principled at every step / structurally-independent reference /
   FP-non-associativity dimensionally explainable). For the
   curvilinear :math:`M_{\rm spatial} = (L+C) - M_{\rm angular\_redist}`
   subtraction, the algebra-decomposition invariant passes at
   ~16×ULP.

3. **Structural-independence ground at L1**. The pre-snapshot
   regression tests are bit-identity against the OLD code; they
   cannot catch a bug that was ALREADY in the old code and survived.
   The L1 ground truth is two-pillared (per
   the ``vv-principles`` skill):

   - **Closed-form pillar**: :math:`k_\infty = \nu\Sigma_f / \Sigma_a`
     on homogeneous reflective slab / sphere / cylinder. Verified at
     ``tests/sn/l1_analytical/test_kinf_homogeneous.py``. This is the
     eigenvalue reference — MMS does NOT prove eigenvalues per
     the ``vv-principles`` skill §"What each pillar
     proves".

   - **MMS pillar**: P1 anisotropic manufactured-source convergence
     at ``tests/sn/test_mms_aniso.py``,
     ``tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py``,
     ``tests/sn/test_mms_heterogeneous.py``, and
     ``tests/sn/test_mms_2d.py``. The MMS source is structurally
     independent of the operator-algebra path (derived by SymPy in
     ``orpheus/derivations/sn/mms_*.py``); it catches flux-shape and
     convergence-order errors that snapshot bit-identity cannot.

4. **Algebraic-identity gates** (new in Wave T). Each touched
   operator passes the algebra contracts:

   - :meth:`TensorProductOperator.assert_separable` passes on every
     TP-shaped operator (BC realizers, fission). This is structurally
     inapplicable to :class:`OperatorSum`-of-bespoke-leaves (T.3
     scattering kernel; the T.4 streaming matvec, fused since #238) —
     see the "out of scope" note in
     ``.claude/plans/wave_t_tensor_network.md`` §6 T.5.

   - **(#238 retired)** the Wave-T algebra-decomposition invariant
     :math:`(L+C)\,\psi \equiv M_{\rm spatial}\,\psi +
     M_{\rm angular\_redist}\,\psi` pinned the typed-leaf split; with the
     split removed the surviving guarantee that the curvilinear
     angular-redistribution term is correct is the anisotropic
     curvilinear MMS (:ref:`sn-mms-curvilinear-aniso-verification`).

   - :math:`(L+C).{\rm solve}(q)` bit-identical pre/post-Wave-T,
     verifying the WDD sweep procedural inverse was NOT touched
     (the :class:`InvertibleOperator.solve` body runs the procedural
     algorithm on the operator's own
     :attr:`~orpheus.sn.operators.streaming.InvertibleOperator.loss_representation`
     since S6.5 — at Wave T it was the free function ``transport_sweep``).

5. **Performance regression gate**. The 1-D slab Krylov benchmark
   measured median 1.04× pre-T.4 baseline (under the 5% threshold).
   (#238 retired the ``M_spatial`` / ``M_angular_redist`` cached
   properties along with the typed-leaf split; the fused matvec carries
   no per-leaf construction cost.)


What :class:`SumOfTensorProductsOperator` was supposed to do — and didn't
-------------------------------------------------------------------------

Grand Report v3 §15.2 (lines 2046-2086) names the canonical scattering
form

.. math::

   S \;=\; \sum_{\ell=0}^{L}
     \Sigma_{s,\ell}\, \otimes\, A_\ell\, \otimes\, G_\ell

with :math:`A_\ell` the angular Pℓ-projection factor,
:math:`\Sigma_{s,\ell}` the per-ℓ group-coupling factor, and
:math:`G_\ell` the per-ℓ spatial factor. Wave T's original T.3 plan
targeted this SOTP form for the scattering kernel.

The design fork (T.3 spec Q6) surfaced that the per-material per-ℓ
einsum in
:meth:`MaterialXSField.apply_legendre_scattering_moments`
**couples the group axis with the spatial axis** — the per-cell
material id ``cells_by_material[mid]`` selects the per-material
scattering matrix :math:`\Sigma_{s,\ell}^{m(\vec r)}`. There is no
factor design where one factor acts on the group axis alone and
broadcasts on the spatial axis; the per-cell material id breaks the
broadcast contract.

The user-resolved math-honest fallback shipped at commit ``9f85c5d``:
:class:`OperatorSum` over per-ℓ :class:`_PerLegendreOrderScattering`
bespoke leaves. The §15.2 *form* is preserved at the summation level
(one summand per Legendre order); the per-summand decomposition into
:math:`R_\ell \circ \Lambda_\ell \circ M_\ell` is a procedural
composition, not a tensor product.

The same master condition applies to T.4-spatial (per-direction WDD
recurrence) and T.4-curvilinear (M-M half-grid recurrence). Two of
the three originally-SOTP-targeted Wave-T substeps fell back to
:class:`OperatorSum`-of-bespoke-leaves; only T.1 (BC realizers) and
T.2 (fission rank-1) cleanly support the TP form.

**Implication for Wave O (Issue #208)**. The operator-role typing
work MUST accommodate non-SOTP :class:`OperatorSum` summands. Any
contract of the form "every BulkOperator summand IS a
:class:`TensorProductOperator`" forecloses scattering, streaming
spatial, and curvilinear angular redistribution. The five-shape
catalogue in :ref:`wave-t-shape-table` is the constraint the Wave O
typing must respect.


Cross-references
----------------

- **Wave T plan** (canonical reference for substep sequencing,
  architectural decisions, deviations from §15.2):
  ``.claude/plans/wave_t_tensor_network.md``.
- **T.4 verification spec** (Q1-Q5 architectural decisions, risk
  register, test catalogue):
  ``.claude/agent-memory/test-architect/wave_t_t4_streaming_verification_spec.md``.
- **Grand Report v3** §15 (V = X ⊗ Ω ⊗ G), §15.1 (streaming as sum of
  tensor products), §15.2 (scattering as sum of tensor products),
  §16A.10 (BC as tensor network), §35 (commandments), north-star line
  5697.
- **Shipped commits**: ``fa13e78`` (T.1 BC), ``0b2848b`` (T.2
  fission), ``9f85c5d`` (T.3b kernel), ``03bcdba`` (T.3c
  build_aniso_source rewire), ``cb18fdb`` (T.4a snapshots),
  ``c55b505`` (T.4b slab M_spatial), ``90e7d4e`` (T.4c curvilinear
  M_angular_redist).
- **Code anchors**:

  - :class:`orpheus.numerics.operator.TensorProductOperator`,
    :class:`orpheus.numerics.operator.SumOfTensorProductsOperator`,
    :class:`orpheus.numerics.operator.OperatorSum`,
    :class:`orpheus.numerics.operator.RankOneOperator`,
    :class:`orpheus.numerics.operator.IdentityOperator`,
    :class:`orpheus.numerics.operator.ZeroOperator`.
  - :class:`orpheus.sn.boundary.realizer.SNBoundaryRealizer` —
    the BC realizer dispatching the T.1 lifts.
  - :class:`orpheus.transport.operators.fission.FissionOperator` and its
    :attr:`~orpheus.transport.operators.fission.FissionOperator.kernel` property.
  - :class:`orpheus.transport.operators.scattering.ScatteringOperator` and its
    :attr:`~orpheus.transport.operators.scattering.ScatteringOperator.kernel` /
    :attr:`~orpheus.transport.operators.scattering.ScatteringOperator.kernel_summands`
    properties; the bespoke
    :class:`orpheus.transport.operators.scattering._PerLegendreOrderScattering` leaf.
  - :class:`orpheus.sn.operators.streaming.StreamingOperator` and its
    :meth:`~orpheus.sn.operators.streaming.StreamingOperator.apply` /
    :meth:`~orpheus.sn.operators.streaming.StreamingOperator.apply_transpose`
    public matvec surface. (#238 retired the per-direction
    ``M_spatial`` / ``M_angular_redist`` typed-leaf split — the
    ``_SpatialSweepDirection`` / ``_MSpatialOperatorSum`` /
    ``AngularRedistributionOperator`` leaves had no production
    consumer.)
  - :meth:`orpheus.sn.loss_representation._OneDimScanWalk._apply_walk`
    — the fused single-emission 1-D matvec body (``(L+C)ψ``), the
    apply-direction twin of the sweep. #206 Phase C moved the walk off
    the operator INTO the representation; #238 removed the dual-emission
    ``(M_spatial, M_angular_redist)`` arm (no production consumer). The
    public surface is :meth:`StreamingOperator.apply`.
  - :class:`orpheus.sn.spatial.pole_angular_closure.PoleAngularClosureBase`
    — the M-M closure data and per-cell algebra primitive.
  - :func:`orpheus.transport.spatial.cell_balance.cell_balance_for_streaming`
    — the three-term cell-balance primitive.


.. _bc-extraction:

Boundary-condition extraction — :math:`B` as a sibling operator (Wave O)
========================================================================

Wave O step O.4a.2 (`Issue #208
<https://github.com/deOliveira-R/ORPHEUS/issues/208>`_, three commits
``d7e1316`` / ``4c0ff96`` / ``2bdc66d``, 2026-06-03) made the realized
boundary law :math:`B` a **first-class sibling** of the streaming +
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
prediction (:ref:`wave-t-tensor-network`): "Wave O typing must accept
non-SOTP summands." :class:`~orpheus.sn.operators.boundary.SNBoundaryOperator`
is the :math:`A_{ss}` leaf — a bespoke
:class:`~orpheus.numerics.operator.LinearOperator` carrying
:attr:`~orpheus.numerics.operator.BlockRole.BOUNDARY`, NOT a
:class:`~orpheus.numerics.operator.TensorProductOperator`.


The block matrix
----------------

On :math:`V = V_{\rm bulk} \oplus V_{\rm boundary}` the four operator
families occupy disjoint blocks of the :math:`2\times 2` block matrix
(:math:`b` = bulk row/column, :math:`s` = surface/trace row/column):

.. math::
   :label: bc-extraction-block-matrix

   \underbrace{
   \begin{bmatrix} A_{bb} & 0 \\ 0 & 0 \end{bmatrix}
   }_{C,\,S,\,F\;(\text{BULK})}
   \;+\;
   \underbrace{
   \begin{bmatrix} A_{bb} & A_{bs} \\ A_{sb} & 0 \end{bmatrix}
   }_{L_{\rm full}\;(\text{FULL})}
   \;-\;
   \underbrace{
   \begin{bmatrix} 0 & 0 \\ 0 & A_{ss} \end{bmatrix}
   }_{B\;(\text{BOUNDARY})}

.. vv-status: bc-extraction-block-matrix documented

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
       via :class:`~orpheus.sn.operators.streaming.InvertibleOperator` ``L+C``)
     - ``FULL``
     - Reads :math:`\psi.\text{bulk}` **and** the *given*
       :math:`\psi.\text{inflow}` trace; writes :math:`\psi.\text{bulk}`
       and the :math:`\psi.\text{outflow}` trace.
     - :math:`A_{bb}` (streaming) + :math:`A_{bs}` (inflow seeds the
       sweep) + :math:`A_{sb}` (sweep produces outflow). The
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
       realized per-face law :math:`\psi.\text{inflow} =
       B\,\psi.\text{outflow}`.
     - :math:`A_{ss}` only; emits on the **inflow row ONLY**
       (see design correction 2 below).

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

For vacuum :math:`B = 0` (no inflow); for reflective/white/albedo
:math:`B` is the realized :math:`R\,G` reflection that the
pre-extraction sweep applied at the boundary edge.


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

**2.** :math:`B` **must project to the inflow row.**
The realized per-face law is a **full-face operator**: a specular
:class:`~orpheus.numerics.operator.PermutationOperator` for reflective,
an :class:`~orpheus.sn.boundary.angular.AngularAverageOperator` for
white. Its permutation maps the input's *inflow* slots onto the
*output's outflow* slots (a spurious :math:`R\cdot\psi.\text{inflow}`),
because the permutation is defined on the whole face, not just the
:math:`A_{ss}` :math:`V_{\rm outflow} \to V_{\rm inflow}` map. In the
legacy sweep this was harmless — the sweep only ever read the
inflow slots of ``bc.apply(face)``, discarding the outflow output.
But as a sibling :math:`-B` on the direct-sum state, a non-zero
outflow emission would corrupt the outflow-definition residual
:math:`r_{\rm outflow}` (which must carry **no** :math:`B` term). The
fix:
:meth:`SNBoundaryOperator._apply_faces <orpheus.sn.operators.boundary.SNBoundaryOperator>`
**projects** the emission onto the codomain row — ``apply`` writes the
``inflow_indices_for_face`` slots; ``apply_transpose`` writes the
``outflow_indices_for_face`` slots. *Empirically confirmed before the
fix*: the outflow slots carried nonzero :math:`R\cdot\psi.\text{inflow}`.

**3. The bare sweep seeds inflow from** :math:`\text{rhs.boundary}`,
**not** :math:`\text{initial\_guess.boundary}`.
Under the extraction the WDD sweep ``(L+C).solve`` is **bare** (see
:ref:`bare-sweep-extraction` in :doc:`discrete_ordinates`): it reads
the seeded inflow trace directly instead of re-applying ``bc``.
:meth:`InvertibleOperator._solve_timed_full_field <orpheus.sn.operators.streaming.InvertibleOperator._solve_timed_full_field>`
must therefore seed the sweep's boundary buffer from
:math:`\text{rhs.boundary}` — the *boundary source*
:math:`q.\text{boundary} + B\,\psi.\text{outflow}` — **not** from the
iterate ``initial_guess.boundary`` (the retired partner-flux carrier).
The ``initial_guess`` still threads the bulk Carlson warm start;
only the boundary seed moved.


.. _bc-extraction-variadic-driver:

The honest :math:`L+C-S-B` driver via variadic couplings (Wave O step O.2a)
---------------------------------------------------------------------------

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

   (L+C,\; S,\; B)
   \quad\Longrightarrow\quad
   \underbrace{(L+C).\text{apply} - S.\text{apply}
               - B.\text{apply}}_{\equiv\,(L+C-S-B)\,\psi}
   \,,\qquad
   \text{rhs} = q_{\rm ext} + S.\text{apply}(\psi) + B.\text{apply}(\psi)

.. vv-status: bc-extraction-within-group-decomposition documented

assembled once by the single-source-of-truth builder
:func:`~orpheus.sn.coupled_system.build_within_group_system`, which returns
the frozen :class:`~orpheus.sn.coupled_system.WithinGroupSystem` record —
the loss grid together with its **named regular splitting**
:math:`A = M - N` (Hackbusch 2016 §11). On a seedless (slab / cylinder /
Cartesian) mesh the record degrades to exactly this triple: its ``resolvent``
is :math:`M = (L+C)` — the invertible resolvent
(:class:`~orpheus.sn.operators.streaming.InvertibleOperator`, ``.solve`` = the WDD
sweep) — and its ``gains`` are :math:`N = (S,\ B_a)`, the two lagged
couplings the driver applies: the bulk scattering gain
(:class:`~orpheus.transport.operators.scattering.ScatteringOperator`,
:attr:`block_role <orpheus.numerics.operator.BlockRole>` ``BULK``) and the
boundary reflection gain
(:class:`~orpheus.sn.operators.boundary.SNBoundaryOperator`,
``block_role`` ``BOUNDARY``). On a *carrying* sphere the record poses the
ψ½ System B as the second row of a coupled :math:`M - N` block grid — the
full 2×2 coupled block operator documented in
:ref:`coupled-block-operator` (its starting-direction physics in
:ref:`sn-282-direct-starting-direction-solve` in :doc:`discrete_ordinates`);
the former ``_within_group_triple`` / ``_lagged_gains`` construction pair
retired into this one builder at B.2d. The :math:`B\,\psi.\text{outflow}`
term lands on :math:`\text{rhs.boundary}`, which the bare ``(L+C).solve``
sweep reads as the inflow seed (:ref:`bare-sweep-extraction` in
:doc:`discrete_ordinates`).

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

#. **The adjoint metric lives on the trace.** :math:`B` lives on the
   boundary trace (:attr:`domain <orpheus.sn.operators.boundary.SNBoundaryOperator.domain>`
   ``= sn_mesh.angular_trace``), and the cosine-weighted
   :math:`|\Omega\cdot\hat n|\,w` adjoint metric (Wave O step O.2 — the
   codomain inner product of :math:`L`'s boundary-trace block) lives on
   that **trace** domain, not the bulk. Folding :math:`B` into the
   bulk :math:`S` would erase the trace typing the future adjoint
   ``.H`` needs.
#. :math:`B` **cannot join the** :math:`L+C` **preconditioner.** A
   generic :class:`~orpheus.numerics.operator.OperatorSum` carries **no**
   direct-sweep ``solve`` — its inverse action is the *iterative*
   :class:`~orpheus.numerics.green_operator.GreenOperator` splitting
   (:ref:`green-operator`), not the ``(L+C)``
   :class:`~orpheus.sn.operators.streaming.InvertibleOperator` subclass's
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
   *honest* :math:`(L+C-S-B)\,\psi` and :math:`q_{\rm ext}+S\psi+B\psi`
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
   * - **Direct helper**
       (:func:`_reflect_outflow_into_inflow <orpheus.sn.solver._reflect_outflow_into_inflow>`)
     - Fills each face's inflow slots with
       :math:`B\,\psi.\text{outflow}` in place on the boundary buffer,
       via the same :class:`SNBoundaryOperator`, before the bare
       sweep.
     - The loops that have **no driver to route through**: the final
       eigenvalue reconstruction sweep in
       :func:`solve_sn <orpheus.sn.solver.solve_sn>`, and the
       octant-restricted Gauss-Seidel variant (Phase 3). The direct
       fixed-source SI loop now routes through the variadic driver, so
       it no longer needs this helper.

The direct helper is **not** a fold of :math:`B` into :math:`S`: it is
the trace-only :math:`A_{ss}` action of the *same* :math:`B`, expressed
on the boundary trace alone. Both routes therefore deliver the
identical :math:`-B` coupling, and cannot drift, because both descend
from :meth:`SNBoundaryOperator._reflect_trace <orpheus.sn.operators.boundary.SNBoundaryOperator>`
(:ref:`bc-extraction-reflect-trace`).


.. _bc-extraction-reflect-trace:

The trace-only :math:`A_{ss}` leaf — :meth:`reflect_into_inflow`
----------------------------------------------------------------

:math:`B` is the :math:`A_{ss}` block :math:`V_{\rm outflow} \to
V_{\rm inflow}`: it maps the *outflow* trace to the *inflow* trace.
Both delivery routes ultimately need the same per-face action — apply
each face's realized law (the specular
:class:`~orpheus.numerics.operator.PermutationOperator` for reflective,
:class:`~orpheus.sn.boundary.angular.AngularAverageOperator` for
white, zero for vacuum) and project onto the inflow row. To guarantee
they cannot drift, that action is the single
:meth:`SNBoundaryOperator._reflect_trace <orpheus.sn.operators.boundary.SNBoundaryOperator>`
core, and both the full-field forward action
:meth:`B.apply <orpheus.sn.operators.boundary.SNBoundaryOperator.apply>`
and the new trace-only leaf
:meth:`B.reflect_into_inflow <orpheus.sn.operators.boundary.SNBoundaryOperator.reflect_into_inflow>`
route through it (Wave O step O.2a, commit ``8563f4b``).

The leaf exists because the direct helper does not need a full field.
:meth:`B.apply` operates on a :class:`~orpheus.transport.full_field.FullField`
(zero bulk, trace populated) — the timeless, history-blind operator carrier
(:meth:`SNBoundaryOperator.apply <orpheus.sn.operators.boundary.SNBoundaryOperator.apply>`
is the base arrow ``FullField -> FullField``; the comonad lives on the
driver), the bulk only a carrier to reach the :math:`A_{ss}` boundary block. The pre-extraction direct helper
fabricated a throwaway zero-bulk field purely to call ``B.apply`` and
then discarded the (zero) bulk output.
:meth:`reflect_into_inflow <orpheus.sn.operators.boundary.SNBoundaryOperator.reflect_into_inflow>`
takes a bare :class:`~orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux`
trace and returns the boundary-only
:class:`~orpheus.transport.source_sinks.angular_boundary_source_sink.AngularBoundarySourceSink`
directly — no zero-bulk probe.

The projection onto the inflow row is load-bearing: the realized law is
a *full-face* operator (the specular permutation also maps the input
inflow slots onto the *outflow* slots, :math:`R\,\psi.\text{inflow}`).
The legacy in-sweep ``bc.apply`` only ever read the inflow slots of its
output, so that spurious outflow emission was harmless. But as the
sibling :math:`-B` reading the *whole* boundary block, a non-zero
outflow emission would corrupt the outflow-definition residual
:math:`r_{\rm outflow}` (which must carry **no** :math:`B` term —
:ref:`bc-extraction-design-corrections`). So
:meth:`_reflect_trace <orpheus.sn.operators.boundary.SNBoundaryOperator>`
projects the forward action onto ``inflow_indices_for_face`` and the
Euclidean transpose onto ``outflow_indices_for_face``.


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
geometry paths cannot drift: ``sn_mesh.reduced is not None`` is the
**same** predicate :func:`~orpheus.sn.loss_representation.transport_sweep` uses to
select the 1-D scan vs the 2-D wavefront body, and the **same**
predicate the direct-helper guards
(:func:`_solve_fixed_source_si <orpheus.sn.solver._solve_fixed_source_si>`,
:func:`solve_sn <orpheus.sn.solver.solve_sn>`) check before calling
:func:`_reflect_outflow_into_inflow <orpheus.sn.solver._reflect_outflow_into_inflow>`.
Both branches are now bare-sweep + sibling :math:`-B`; the predicate
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
:ref:`bc-extraction-g-adjoint`). A consequence is that the original
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

Because the variadic :math:`-S - B` gains ride the **bare** sweep for
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
(:meth:`AngularSourceSink.from_isotropic <orpheus.transport.source_sinks.angular_source_sink.AngularSourceSink.from_isotropic>`
bulk + :meth:`AngularBoundarySourceSink.zeros_on <orpheus.transport.fields._bases.AngularBoundaryField.zeros_on>`
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
loop) by ``tests/sn/eigenvalue/test_keff_2d.py::TestSIKrylov2DEquivalence``:

.. list-table:: 2-D Cartesian SI/Krylov verification (Wave O step #208)
   :header-rows: 1
   :widths: 38 34 28

   * - Leg
     - Reference (pillar)
     - Result
   * - Default-entry homogeneous (1G / 2G / 4G)
       (:func:`test_default_entry_hits_kinf <tests.sn.eigenvalue.test_keff_2d.TestSIKrylov2DEquivalence.test_default_entry_hits_kinf>`)
     - Closed-form :math:`k_\infty = \lambda_{\max}(A^{-1}F)`
       (1g → 1.5, 2g → 1.875, 4g → 1.4878)
     - SI hits :math:`k_\infty` to :math:`< 10^{-8}`
   * - Heterogeneous 2G fuel\|moderator, **non-flat** flux
       (:func:`test_si_krylov_heterogeneous_2g_nonflat_flux <tests.sn.eigenvalue.test_keff_2d.TestSIKrylov2DEquivalence.test_si_krylov_heterogeneous_2g_nonflat_flux>`)
     - SI vs Krylov flux **shape** + eigenvalue
     - flux shape agrees to :math:`\sim 10^{-9}`
   * - 2-D SI :math:`k_{\rm eff}` Cauchy under refinement
       (:func:`test_si_2d_keff_converges_under_refinement <tests.sn.eigenvalue.test_keff_2d.TestSIKrylov2DEquivalence.test_si_2d_keff_converges_under_refinement>`)
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
(``tests/sn/sweep/cartesian_2d/test_discrete_ordinates_2d.py::test_do_mesh_convergence``,
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
vacuum path is **bit-identical** to the pre-extraction matvec. Verified
by:

- the matvec 18-baseline snapshot
  (:func:`np.array_equal` against the pre-O.4a.2 captures across
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
       (``tests/sn/spatial/test_streaming_equilibrium_curvilinear.py``)
     - Analytical infinite-medium balance
       :math:`\phi = q/\Sigma_t` (closed-form)
     - ``source_iteration`` AND ``krylov``
   * - Reflective :math:`k_\infty` homogeneous
       (``tests/sn/l1_analytical/test_kinf_homogeneous.py``)
     - :math:`k_\infty = \nu\Sigma_f / \Sigma_a` (closed-form
       eigenvalue — MMS does NOT prove eigenvalues)
     - both
   * - ``test_si_carve_recovers_analytical_kinf``
       (``tests/sn/operators/test_invertible_operator.py``)
     - Analytical :math:`k_\infty` via the SI path with :math:`B`
       folded (closed-form)
     - SI path
   * - Invertible-operator :math:`Q/\Sigma_t` recovery
       (``tests/sn/operators/test_invertible_operator.py``)
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
retired fold) to :math:`(L+C) - S - B` (the two separate gains of
:ref:`bc-extraction-variadic-driver`), and the rhs symmetrically,
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

Operator-output role typing — :math:`A\psi` is a source/sink (Wave O step B.5.2)
--------------------------------------------------------------------------------

Wave O step B.5.2 (`Issue #208
<https://github.com/deOliveira-R/ORPHEUS/issues/208>`_, commit
``6ef5063``, 2026-06-03) retyped every SN operator's ``.apply`` output
``.boundary`` from
:class:`~orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux` to
:class:`~orpheus.transport.source_sinks.angular_boundary_source_sink.AngularBoundarySourceSink`
— the *source/sink* role leaf. This completes the **boundary** half of
the B.5 operator-output "dimensional-sin" carve; the **bulk** half
(``.apply.bulk`` →
:class:`~orpheus.transport.source_sinks.angular_source_sink.AngularSourceSink`)
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
     - :class:`~orpheus.transport.source_sinks.angular_source_sink.AngularSourceSink`
       (``f400743``)
     - :class:`~orpheus.transport.fields.angular_flux.AngularFlux`
     - :class:`~orpheus.transport.residuals.angular_residual.AngularResidual`
       (consumed by :func:`~orpheus.sn.solver.evaluate_residual`,
       O.2 — :ref:`affine-typed-residual`)
   * - **boundary** (:math:`V_{\rm inflow} \oplus V_{\rm outflow}`)
     - :class:`~orpheus.transport.source_sinks.angular_boundary_source_sink.AngularBoundarySourceSink`
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
consumer is the honest :math:`L+C-S-F-B` driver of Wave O step
**O.2** (see :ref:`bc-extraction-operator-output-o2`). That consumer
has since landed:
:func:`~orpheus.sn.solver.evaluate_residual` types the balance defect
:math:`(L+C-S-B)\psi - q` via ``from_balance`` (Wave O step O.2
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
   satisfy.** The realized boundary law
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
        - :math:`(L+C).\text{apply} - S.\text{apply} - B.\text{apply}`
        - a **residual** term (subtracted from the diagonal)
      * - SI rhs
        - :math:`q_{\rm ext} + S.\text{apply} + B.\text{apply}`
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
:meth:`InvertibleOperator._solve_timed_full_field <orpheus.sn.operators.streaming.InvertibleOperator._solve_timed_full_field>`,
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

The change is **type-only**: :meth:`AngularBoundarySourceSink.zeros_on <orpheus.transport.fields._bases.AngularBoundaryField.zeros_on>`
and the per-face-view writes produce **bit-identical** ``.values`` —
only the wrapping role-type differs. The dead :class:`AngularBoundaryFlux`
runtime imports were retired from the retyped sites.


.. _bc-extraction-operator-output-o2:

What remains for Wave O step O.2
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Wave O step **O.2a has landed the honest** :math:`L+C-S-B` **driver**
via the variadic couplings of :ref:`bc-extraction-variadic-driver`:
the transitional :math:`S + B` fold is **retired** and :math:`B` is now
a first-class coupling gain. Of the items B.5.2 left for the rest of
O.2, **the adjoint metric and its gate landed in step O.2b
R5** (:ref:`bc-extraction-g-adjoint`), and **the residual column landed
in the O.2 close-out** (:ref:`affine-typed-residual`):

* :meth:`AngularResidual.from_balance <orpheus.transport.residuals.angular_residual.AngularResidual.from_balance>`
  and
  :meth:`AngularBoundaryResidual.from_balance <orpheus.transport.residuals.angular_boundary_residual.AngularBoundaryResidual.from_balance>`
  now have their first production-reachable consumer:
  :func:`~orpheus.sn.solver.evaluate_residual` types the within-group
  balance defect :math:`r = (L+C-S-B)\psi - q` as the composite
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
  V_{\rm trace}` (:ref:`bc-extraction-g-adjoint`). This is exactly why
  :math:`B` stays trace-typed as a separate gain
  (:ref:`bc-extraction-variadic-driver`): a bulk-folded :math:`B` would
  erase the trace metric the adjoint needs.
* **Gate-1.3** (the O.2 adjoint verification gate) **landed in R5** —
  the dense-probe oracle + the L11 wrong-metric control
  (:ref:`bc-extraction-g-adjoint`).

The direct-helper
:func:`_reflect_outflow_into_inflow <orpheus.sn.solver._reflect_outflow_into_inflow>`
also survives O.2a (the driver no longer routes through it, but the
final eigenvalue reconstruction sweep and the Gauss-Seidel variant
still do — :ref:`bc-extraction-two-routes`); the
:attr:`ScatteringOperator.domain` typing completion that was once a
documented seam **landed in P4.5 W-D** — :math:`S` now carries the
composite full-field space (:ref:`bc-extraction-scope-future`).

The residual column is now wired: :class:`AngularBoundaryResidual` and
:class:`~orpheus.transport.residuals.angular_residual.AngularResidual`
are consumed by :func:`~orpheus.sn.solver.evaluate_residual` (Wave O
step O.2 close-out — :ref:`affine-typed-residual`), which also
completes the affine flux algebra (the iterate increment is now the
typed :class:`~orpheus.transport.displacements.angular_displacement.AngularDisplacement`,
and ``flux + flux`` is a :class:`TypeError`). That section is the
canonical home for the **affine-typed triad** (state / displacement /
residual).


.. vv-status note: the operator-output role typing of B.5.2 is a
   type-only refactor (bit-identical ``.values``); its correctness is
   verified by the same gates that verify the extraction
   (:ref:`bc-extraction-numerical-evidence`) plus the type-residual
   gates catalogued below. B.5.2's verification ground:

* the **B0 decision instrument** (``diag_b52_boundary_typing_decision.py``)
  proving OPT-BSS closes both sums while OPT-BR throws the two-hat;
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


.. _affine-typed-field-algebra:

The affine-typed SN field algebra — state, displacement, residual (Wave O step O.2 close-out)
=============================================================================================

Wave O step **O.2** (`Issue #208
<https://github.com/deOliveira-R/ORPHEUS/issues/208>`_, which subsumes
`Issue #201 <https://github.com/deOliveira-R/ORPHEUS/issues/201>`_,
2026-06-08) closes the operator-algebra-completeness wave by making the
SN **flux algebra** as type-honest as the operator algebra above it.
The field-role grid is completed with a fourth column — the
**displacement** — and the residual column gets its first production
consumer.

.. admonition:: Key Facts (the affine triad)
   :class: tip

   - **Three quantities, three roles.** The flux **state** :math:`\psi`
     (:math:`[1/(\mathrm{cm^2 \cdot s \cdot sr})]`), the flux
     **displacement** :math:`\Delta\psi` (same units / shape / metric,
     a DISTINCT role: a *step*, not a *state*), and the equation
     **residual** :math:`r = A\psi - q`
     (:math:`[1/(\mathrm{cm^3 \cdot s \cdot sr})]`, a rate-density
     defect).
   - **Two dimensional universes connected by** the loss operator
     :math:`A = L + C`. Flux and rate-density are joined by :math:`A`
     (:math:`A`: apply, flux → rate-density;
     :math:`A^{-1}`: solve, rate-density → flux — the sweep). State and
     displacement live in the **flux** universe; source and residual
     live in the **rate-density** universe.
   - **Flux states are an affine space** :math:`A` **over a difference
     vector space** :math:`V`. The increment
     :math:`\Delta\psi \in V`; ``flux ⊖ flux → FluxDisplacement`` is the
     only mint; ``flux ⊕ displacement → flux`` is the torsor action;
     ``flux + flux`` is a :class:`TypeError` (no origin); only
     ``affine_combination`` (:math:`\sum\lambda_i = 1`) blends fluxes.
   - **The displacement carries the convergence diagnostics a state
     cannot** (``contraction_ratio`` :math:`\approx c`,
     ``true_error_estimate`` = the :math:`c\to 1` false-convergence
     fix, ``where_largest`` = the per-cell convergence map).
   - **The typed residual** is minted by the named composition
     :meth:`~orpheus.transport.residuals.angular_residual.AngularResidual.from_balance`
     and consumed by :func:`~orpheus.sn.solver.evaluate_residual` (the
     #208 box-7 consumer) — an additive diagnostic + the consistent-DSA
     (Issue #2) substrate, NOT in the convergence path.
   - **Zero numerical change.** The whole carve is type-only on the
     converged answer: the SI stopping norm stays the flat
     :math:`\lVert\Delta\psi\rVert`, the residual evaluation is
     additive — so the converged flux is **bit-identical** to the
     pre-carve reference (frozen ``sha256``).


The two dimensional universes — flux and rate-density
-----------------------------------------------------

Every typed SN field is a *quantity* with a *dimension* (the View-G
decision, Issues #205 / #207: units live on the field, not the
:class:`~orpheus.numerics.space.FunctionSpace`; see
:class:`~orpheus.numerics.field.Field`). The transport solve moves
between exactly two dimensional universes, connected by the loss
operator :math:`A = L + C`:

.. list-table:: The two universes connected by :math:`A` / :math:`A^{-1}`
   :header-rows: 1
   :widths: 22 26 26 26

   * - Universe
     - Units
     - State role
     - Defect / source role
   * - **Flux**
     - :math:`1/(\mathrm{cm^2 \cdot s \cdot sr})`
     - :class:`~orpheus.transport.fields.angular_flux.AngularFlux`
       :math:`\psi`, and its tangent
       :class:`~orpheus.transport.displacements.angular_displacement.AngularDisplacement`
       :math:`\Delta\psi`
     - —
   * - **Rate-density**
     - :math:`1/(\mathrm{cm^3 \cdot s \cdot sr})`
     - —
     - :class:`~orpheus.transport.source_sinks.angular_source_sink.AngularSourceSink`
       :math:`q`, :math:`A\psi`; and
       :class:`~orpheus.transport.residuals.angular_residual.AngularResidual`
       :math:`r`

The map between them is dimensional: applying :math:`A` (which carries
:math:`\hat\Omega\cdot\nabla + \Sigma_t`, both
:math:`1/\mathrm{cm}`) sends a flux to a rate-density (``A.apply`` →
:class:`AngularSourceSink`, documented at
:ref:`bc-extraction-operator-output-typing`); inverting it
(``A.solve``, the sweep :math:`A^{-1} = (L+C)^{-1}`) sends a rate-density
back to a flux (``.solve`` → :class:`AngularFlux`). The source side (the B.5.2
operator-output typing) and the residual side are BOTH rate-density;
the state side and the displacement side are BOTH flux. This is the
load-bearing distinction the role grid encodes: *same universe* grants
permission to add in linear algebra, but the *class* still gates
meaning (a residual and a source share rate-density units yet are
different classes — see :class:`~orpheus.numerics.field.Field`).


The affine structure of flux states
------------------------------------

The new core of step O.2 is the recognition that the flux **states**
do not form a vector space — they form an **affine space**. An affine
space :math:`A` has points but no distinguished origin; differences of
points live in an associated **difference vector space** :math:`V`
(the *torsor* :math:`A` over :math:`V`). The torsor axioms ARE the
field algebra (the SymPy-free derivation lives in
:class:`~orpheus.transport.fields._flux_role.FluxRole`):

.. math::
   :label: affine-torsor-algebra

   \psi_2 \ominus \psi_1 &\;\to\; \Delta\psi \in V
       \qquad&&\text{(the iterate increment; the ONLY mint of a displacement)}\\
   \psi \oplus \Delta\psi &\;\to\; \psi' \in A
       \qquad&&\text{(the torsor action } A\times V\to A\text{; the update step)}\\
   \psi_1 \oplus \psi_2 &\;\to\; \bot
       \qquad&&\text{(no origin — } \textstyle\sum\lambda_i = 2 \text{ lands off } A)\\
   \textstyle\sum_i \lambda_i \psi_i,\ \textstyle\sum_i\lambda_i=1
       &\;\to\; \psi \in A
       \qquad&&\text{(affine combination — the only legal multi-flux blend)}

.. vv-status: affine-torsor-algebra documented

Reading the four lines:

- **The mint** ``flux ⊖ flux → FluxDisplacement``
  (:meth:`FluxRole.__sub__ <orpheus.transport.fields._flux_role.FluxRole.__sub__>`).
  The difference of two flux states is NOT a state — it is the vector
  between two points, an element of :math:`V`. It is the iterate
  increment :math:`\Delta\psi`. This is the *only* way a displacement
  is born (mirroring the residual, which is born *only* from
  ``from_balance``).

- **The torsor action** ``flux ⊕ displacement → flux``
  (:meth:`FluxRole.__add__ <orpheus.transport.fields._flux_role.FluxRole.__add__>`).
  A point plus a vector is a point: the relaxation/update step
  :math:`\psi^{(i+1)} = \psi^{(i)} \oplus \Delta\psi`.

- **The forbidden** ``flux + flux → TypeError``. Two points cannot be
  added: an affine combination with :math:`\sum\lambda_i = 2` lands off
  the affine subspace :math:`A`. The :class:`TypeError` names the legal
  alternatives in its message (``affine_combination`` for a blend,
  ``ψ + Δψ`` for an update, or add the rate-density :math:`L`-actions).

- **The affine combination**
  (:meth:`FluxRole.affine_combination <orpheus.transport.fields._flux_role.FluxRole.affine_combination>`).
  :math:`\sum_i\lambda_i\psi_i` is a flux **iff**
  :math:`\sum_i\lambda_i = 1` (the partition-of-unity / affine
  constraint, enforced to :math:`10^{-12}`). The relaxation
  :math:`\omega\,\psi_{\rm new} + (1-\omega)\,\psi_{\rm old}` is the
  canonical instance — its coefficients sum to 1, equivalently the
  torsor form :math:`\psi_{\rm old} \oplus \omega(\psi_{\rm new}\ominus
  \psi_{\rm old})`.

**Scalar scaling is untouched.** ``__mul__`` / ``__rmul__`` /
``__truediv__`` / ``__neg__`` take a SCALAR and are inherited from
:class:`~orpheus.numerics.field.Field`, so eigenvalue normalisation
:math:`\psi/k`, :math:`\omega\cdot\Delta\psi`, and unary negation all
survive. The gate forbids ONLY the binary ``flux + flux``; ``__sub__``
is retyped, and the cross-class / cross-mesh guards are inherited
unchanged.

.. note:: **The torsor round-trip is exact up to FP rounding, NOT
   bit-exact.** The identity :math:`\psi_1 \oplus (\psi_2 \ominus
   \psi_1) = \psi_2` holds in real arithmetic, but in IEEE-754
   :math:`a + (b - a) \neq b` bit-for-bit (the subtraction rounds), so
   the foundation test
   ``tests/numerics/test_affine_flux_algebra.py::test_mint_stores_raw_difference_and_torsor_round_trips``
   asserts the round-trip to ``nulp=8`` (8 ULP), while the *individual*
   mint (``Δψ.values == ψ₂.values − ψ₁.values``) and the *individual*
   torsor add (``(ψ₁⊕Δψ).values == ψ₁.values + Δψ.values``) are
   bit-exact. This is the bit-identity-vs-principled-equivalence
   distinction applied at the algebra level.


Why affine — the design rationale
----------------------------------

The pre-carve algebra typed ``AngularFlux + AngularFlux →
AngularFlux``: it treated :math:`\psi=0` as an additive origin. That is
a literal **affine-axiom violation** — flux states have no natural
zero (the zero flux is a *chosen base point*, not an identity), so
adding two of them is meaningless (it lands off :math:`A`). The same
pre-carve ``__sub__`` mis-typed the increment :math:`\Delta\psi` as a
*state*. This is the **same sin** that Issue #201 fixes for the
residual: the bare subtraction :math:`A\psi - q` would typecheck as
``AngularSourceSink − AngularSourceSink`` yet MIS-TYPE the defect as a
source instead of a residual (see
:ref:`bc-extraction-operator-output-typing`). In both cases a role
transition was hidden inside a same-shape ``-`` that the type system
could not see.

Three expert frames converged on the affine fix (the investigation
history; the structural decision lives here, the full record in
`Issue #208 <https://github.com/deOliveira-R/ORPHEUS/issues/208>`_):

.. list-table:: The three frames that converged on the displacement type
   :header-rows: 1
   :widths: 26 40 34

   * - Frame
     - The structure it sees
     - What it contributes
   * - **Affine geometry / torsor**
     - flux states are points in :math:`A`; :math:`\Delta\psi` is a
       vector in :math:`V`; the three ops are the torsor axioms
       :eq:`affine-torsor-algebra`
     - ``state + state`` becomes UNREPRESENTABLE by construction; the
       #201 gate is a type consequence, not a runtime check
   * - **Banach fixed-point / contraction**
     - :math:`\Delta\psi^{(i+1)} = M\,\Delta\psi^{(i)}`,
       :math:`M = A^{-1}(S+B)` — an exact linear recurrence on the
       increments
     - the displacement is the natural carrier of the contraction
       factor :math:`\rho`, the a-posteriori bound, the Aitken
       extrapolation — data a state has nowhere to put
   * - **Krylov / residual dual**
     - :math:`\Delta\psi = A^{-1}\,\tilde r` — the increment and the
       residual are the SAME defect viewed in the two universes
       connected by :math:`A`
     - the displacement (flux, primal, SI-native) and the residual
       (rate-density, dual, Krylov/DSA-native) are duals, not
       competitors; choosing a stopping criterion is a *frame* choice

The decisive point all three share: **the displacement is the only
object that knows "previous"/"step"**, so it is the natural — and
type-safe — home for the convergence diagnostics that a flux STATE
structurally cannot carry. A flux is a snapshot; it has no notion of an
increment, a contraction factor, or a previous iterate. Typing the
increment as a state both (a) admits the illegal ``state + state`` and
(b) strands the contraction data with no home — exactly what the
pre-carve SI loop did (it computed :math:`\psi - \psi_{\rm prev}`,
took a scalar norm, and **discarded** the increment).

.. note:: **Why mint the full displacement type at one consumer.** The
   project's "wait for ≥2 consumers" rule guards against *premature
   abstraction with no clear benefit*; it does NOT block a genuine
   primitive for expressiveness when the benefit is established. Here
   both expert frames agree the benefit is real (illegal states made
   unrepresentable + the diagnostics get a type-safe home), so the
   displacement type is built fully even though the SI loop is its
   first consumer. See ``feedback_unify_after_two_instances`` (the
   nuance added this session) — this is a deliberate user-steered
   exception to the defer-abstraction default.


The displacement leaves and the role-grid completion
-----------------------------------------------------

The affine algebra is defined ONCE on the
:class:`~orpheus.transport.fields._flux_role.FluxRole` mixin and mixed
into the four flux-STATE leaves
(:class:`~orpheus.transport.fields.angular_flux.AngularFlux`,
:class:`~orpheus.transport.fields.scalar_flux.ScalarFlux`,
:class:`~orpheus.transport.fields.harmonic_moment_flux.HarmonicMomentFlux`,
:class:`~orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux`) — so
every flux leaf gets a displacement sibling for free, mirroring how the
operator block-roles are *derived* not re-declared. The mint copies the
flux's dataclass init-fields (``values`` / ``space`` / ``mesh``, plus
``L`` for the moment family) into the sibling
:class:`~orpheus.transport.displacements._displacement.Displacement`
leaf, swapping ``values`` — uniform across families, and handling the
moment leaf's ``L`` that a ``cls.from_mesh(...)`` engine could not.

.. note:: The interior-face cochain :math:`C^1_{\rm int}`
   (:ref:`wavefront-flux-cochain`) is deliberately **NOT** a
   :class:`FluxRole` leaf — it is a sweep-internal cochain, not an
   iterate state, so even while the ``WavefrontFlux`` type lived it
   stayed a plain :class:`~orpheus.numerics.field.Field` (no
   ``wavefront ± wavefront`` affine algebra). After the carrier retired
   at S6.4(f) the cochain is just plain numpy arrays
   (``_MovingFrontier`` / ``_octant_face_cochain``), so the
   "not a role leaf" conclusion holds *a fortiori* — there is no field
   type at all to mis-place in the role grid.

This makes the field-role grid a clean 4 × 4 (the displacement column
is the affine completion):

.. list-table:: The completed field-role grid (Issue #205 / #201 / #208)
   :header-rows: 1
   :widths: 16 21 21 21 21

   * - Block
     - Flux (state)
     - Source/Sink
     - Residual (defect)
     - Displacement (tangent)
   * - **Angular**
     - :class:`~orpheus.transport.fields.angular_flux.AngularFlux`
     - :class:`~orpheus.transport.source_sinks.angular_source_sink.AngularSourceSink`
     - :class:`~orpheus.transport.residuals.angular_residual.AngularResidual`
     - :class:`~orpheus.transport.displacements.angular_displacement.AngularDisplacement`
   * - **Scalar**
     - :class:`~orpheus.transport.fields.scalar_flux.ScalarFlux`
     - ``ScalarSourceSink``
     - ``ScalarResidual``
     - :class:`~orpheus.transport.displacements.scalar_displacement.ScalarDisplacement`
   * - **Moment**
     - :class:`~orpheus.transport.fields.harmonic_moment_flux.HarmonicMomentFlux`
     - —
     - —
     - :class:`~orpheus.transport.displacements.moment_displacement.MomentDisplacement`
   * - **Boundary**
     - :class:`~orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux`
     - :class:`~orpheus.transport.source_sinks.angular_boundary_source_sink.AngularBoundarySourceSink`
     - :class:`~orpheus.transport.residuals.angular_boundary_residual.AngularBoundaryResidual`
     - :class:`~orpheus.transport.displacements.angular_boundary_displacement.AngularBoundaryDisplacement`

The **displacement column is the dual of the residual column**: a
residual is ``source ⊖ source`` crossing INTO rate-density (a role
*and* a universe change); a displacement is ``flux ⊖ flux`` staying in
flux units but changing role from *state* to *tangent* (a role change
only). ``transport/displacements/`` mirrors ``transport/residuals/``.


The convergence diagnostics the displacement carries
----------------------------------------------------

The diagnostics live on the
:class:`~orpheus.transport.displacements._displacement.Displacement`
mixin (shared by all four displacement leaves). They are the payoff of
the Banach frame: the SI iterate increment obeys
:math:`\Delta\psi^{(i+1)} = M\,\Delta\psi^{(i)}` with
:math:`M = A^{-1}(S+B)`, a contraction with spectral radius
:math:`\rho(M) \le \max_g \Sigma_{s,g}/\Sigma_{t,g} = c`
([AdamsLarsen2002]_).

**The contraction ratio.**
:meth:`~orpheus.transport.displacements._displacement.Displacement.contraction_ratio`
returns the Banach factor

.. math::
   :label: affine-contraction-ratio

   \rho^{(i)} \;\approx\;
   \frac{\lVert \Delta\psi^{(i)} \rVert}{\lVert \Delta\psi^{(i-1)} \rVert} .

.. vv-status: affine-contraction-ratio documented

A value :math:`\rho > 1` diverges (wrong fixed point / unstable
scheme); :math:`\rho \approx 1` is stalled (the :math:`c\to 1` slow
mode; curvilinear / reflective slow modes); :math:`\rho < 1` is
healthy. It turns the :math:`\rho`-blind
:math:`\lVert\Delta\psi\rVert` stopping test honest.

**The true-error estimate — the** :math:`c\to 1` **false-convergence
fix.** This is the load-bearing numerical content. For a geometric
contraction the distance from the current iterate to the fixed point is
the tail sum

.. math::
   :label: affine-true-error

   e^{(i)} \;=\; \sum_{j\ge 0} \rho^{j}\,\lVert\Delta\psi^{(i)}\rVert
          \;=\; \frac{\lVert\Delta\psi^{(i)}\rVert}{1-\rho} ,

.. vv-status: affine-true-error documented

so the bare increment :math:`\lVert\Delta\psi\rVert` **understates** the
true error by a factor :math:`1/(1-\rho)`. As the scattering ratio
:math:`c \to 1` (optically thick, near-pure-scatter), :math:`\rho \to
1` and the understatement blows up: at :math:`c = 0.99`,
:math:`\rho\approx 0.99` and the true error is :math:`\sim 100\times`
the increment — so a solve that "converges" at
:math:`\lVert\Delta\psi\rVert < \text{tol}` is actually
:math:`\sim 100\cdot\text{tol}` from the solution. This is the
canonical source-iteration stall-masking-as-convergence trap
([AdamsLarsen2002]_).
:meth:`~orpheus.transport.displacements._displacement.Displacement.true_error_estimate`
surfaces it (it raises if :math:`\rho \notin [0,1)` — a
non-contracting iteration has no finite geometric-tail estimate). The
displacement is the natural home for this because the estimate needs
both the increment AND :math:`\rho` (which needs the *previous*
increment) — data a flux state cannot hold.

**The convergence map.**
:meth:`~orpheus.transport.displacements._displacement.Displacement.where_largest`
returns the :math:`k` index tuples with the largest
:math:`|\Delta\psi|` — the per-cell / per-group / per-ordinate
convergence map (which entries are not converging: a pole-cell
resonance, a material-interface slow mode, a lagging group).

**Exposed by the SI loop as a debug hook, not the stopping test.**
:class:`~orpheus.numerics.iteration.SourceIteration` records
``self.contraction_ratios`` (the :math:`\rho` history) and
``self.last_displacement`` (the typed leaf, for ``where_largest`` /
``true_error_estimate``) — O(1) extra field memory (one retained
previous displacement). The stopping test is **unchanged**: it stays
the flat relative residual
:math:`\lVert\Delta\psi\rVert / \lVert\psi\rVert` computed with the
flat ravelled Euclidean norm
:meth:`~orpheus.numerics.iteration._l2_norm`, NOT the space-induced
metric :attr:`Field.l2 <orpheus.numerics.field.Field.l2>`. The metric
norm is a **diagnostic only** — switching the stopping norm would
re-interpret ``conv_tol`` and break bit-identity. The synthetic L0
case (a bare ndarray increment) gets no diagnostics and is unaffected.


.. _affine-typed-residual:

The typed equation residual — the box-7 ``from_balance`` consumer
-----------------------------------------------------------------

The residual column of the role grid (B.5.2 left it
*minted-but-consumerless*) is now wired by
:func:`~orpheus.sn.solver.evaluate_residual`, the **#208 box-7
consumer**. It evaluates the within-group balance defect

.. math::
   :label: affine-typed-residual-eq

   r \;=\; (L + C - S - B)\,\psi \;-\; q

.. vv-status: affine-typed-residual-eq documented

as the typed composite
``FullField(bulk=AngularResidual, boundary=AngularBoundaryResidual)`` — the
**timeless** carrier (a residual is a one-shot balance defect, history-free;
P4.5 W-C confines the timed type to the driver iterate) —
minted via the named composition
:meth:`AngularResidual.from_balance <orpheus.transport.residuals.angular_residual.AngularResidual.from_balance>`
/
:meth:`AngularBoundaryResidual.from_balance <orpheus.transport.residuals.angular_boundary_residual.AngularBoundaryResidual.from_balance>`
— **NOT** a bare cross-class ``−`` (which would mis-type the defect as a
source). The operator output :math:`(L+C-S-B)\psi` is a source/sink
composite (the B.5.2 typing); subtracting the source :math:`q` is a
role transition (source operands → residual result), so it must go
through ``from_balance``, exactly as the affine increment must go
through ``__sub__``.

The residual carries three diagnostics (a typed defect a flux cannot
be):

- :meth:`AngularResidual.balance_map <orpheus.transport.residuals.angular_residual.AngularResidual.balance_map>`
  — the per-cell / per-group transport-balance violation
  :math:`\max_n |r_n(\vec r, g)|`. This is the **typed form of the
  per-ordinate flat-flux residual probe** (vv-principles Signature 1 /
  §H3): it exposes WHERE the discrete balance is violated per cell, the
  localised defect that *global* conservation (telescoping particle
  balance) HIDES — for example the curvilinear pole-cell spike (the
  ERR-026 failure-mode-7 class, where the SI gives a large pole-cell
  error while global balance looks fine).
- :func:`~orpheus.sn.solver.boundary_vs_interior_split` — splits the
  composite into :math:`(\lVert r_\partial \rVert, \lVert r_{\rm bulk}
  \rVert)` (with :math:`\sqrt{b^2 + i^2} = \lVert r\rVert`, the same
  flat metric the SI test uses), discriminating a BC-realizer /
  reflective-trace defect from an interior-streaming defect — free from
  the typed composite.
- :meth:`AngularResidual.relative_to <orpheus.transport.residuals.angular_residual.AngularResidual.relative_to>`
  — the tolerance-portable relative residual :math:`\lVert r\rVert /
  \lVert q\rVert` (the bare residual has rate-density units, so its
  magnitude scales with :math:`\Sigma_t \cdot V`; dividing by the drive
  makes it problem-portable).

The residual is **additive / diagnostic** — never in the convergence
path. The SI stopping test stays :math:`\lVert\Delta\psi\rVert`; the
GMRES defect stays the *flat* :math:`b - A\psi` on the raveled vector
(never typed as a field). The typed residual is also the literal
substrate the **consistent DSA** low-order correction (Issue #2) will
consume — DSA computes the transport residual, then a diffusion
correction through the in-algebra diffusion operator
:math:`A_{\rm diff} = L + C - S - B` (now built, #290 P4;
:doc:`/theory/diffusion_1d`); :eq:`affine-typed-residual-eq` is that
transport residual, typed. (``as_dsa_source`` lands WITH DSA #2, per
the build-the-genuine-primitive-defer-the-speculative-tail discipline.)


Numerical evidence
------------------

The carve is verified by three independent gates (all in the worktree,
run under ``-O``):

.. list-table:: Affine-algebra verification gates
   :header-rows: 1
   :widths: 30 22 48

   * - Gate (test)
     - Level / pillar
     - What it proves
   * - ``test_affine_carve_bit_identity.py`` (``foundation``)
     - regression by inheritance
     - The converged ``psi`` / ``phi`` bytes hash to a frozen
       ``sha256`` taken at the pre-carve commit ``63719a2``, for the
       windowed-moment SI path (``si_2d_p1_aniso_het``), the
       full-angular Krylov path (``krylov_2d_p1_aniso_het``), and the
       1-D ``AngularFlux``-bulk SI path (``si_slab_2g_het``). A drift
       here means the carve changed the numerics — sub-ULP, sharper
       than the ``≈1e-11`` DD regression snapshots (which already
       pre-drift ~6920 ULP from Phase 5b/5c).
   * - ``test_flux_displacement_diagnostics.py`` (``l1``)
     - closed-form (:math:`\rho = c`)
     - ``SourceIteration.contraction_ratios`` :math:`\to \rho \approx
       c` on a homogeneous slab: measured :math:`\rho \in [0.40, 0.56]`
       at :math:`c = 0.5` and :math:`\rho \in [0.80, 0.92]` at
       :math:`c = 0.9`, with :math:`\rho(0.9) > \rho(0.5) + 0.2` (the
       discriminating cross-check that :math:`\rho` *tracks* :math:`c`,
       not a constant). 1-group is acceptable HERE because
       :math:`\rho = c` is a flux-shape-independent *rate* claim, NOT
       an eigenvalue / flux-shape claim.
   * - ``test_typed_residual_evaluation.py`` (``l0`` + ``foundation``)
     - term verification
     - ``balance_map`` :math:`\approx 0` at the SI fixed point
       (relative defect :math:`< 10^{-7}`), detectably :math:`\neq 0`
       and localised when one interior cell is perturbed by 10 %
       (:math:`> 100\times` the converged peak, at the perturbed
       cell); ``from_balance`` mints the correct type / units / space
       (and raises on a flux operand);
       :math:`\sqrt{b^2+i^2} = \lVert r\rVert`;
       ``relative_to(q) = \lVert r\rVert/\lVert q\rVert``.

The torsor algebra itself (the mint, the gate, the affine combination,
the telescoping :math:`(\psi_2\ominus\psi_1)+(\psi_3\ominus\psi_2) =
\psi_3\ominus\psi_1`, the round-trip to 8 ULP) is pinned across all
four flux leaves by ``tests/numerics/test_affine_flux_algebra.py``
(``foundation`` — the structural ground is numpy ``+`` / ``-`` on the
same buffer; every negative test is L11-paired with a positive one).


.. _bc-extraction-g-adjoint:

The composite metric-correct G-adjoint — ``op.H`` over ``FullFieldSpace`` (Wave O step O.2b R5)
===============================================================================================

Wave O step **O.2b R5** (`Issue #208
<https://github.com/deOliveira-R/ORPHEUS/issues/208>`_, commits
``89b2f62`` / ``0efd233`` / ``5c06196``, 2026-06-05) discharges the
second open item of :ref:`bc-extraction-operator-output-o2`: it makes
``op.H`` — the Hilbert adjoint of an SN operator composite — the
**metric-correct G-adjoint**

.. (vv-status rationale) The G-adjoint defining identity. The
   verifiable claim is the reciprocity ⟨Aψ,φ⟩_G = ⟨ψ,A†φ⟩_G plus the
   block-diagonal G-fold op.H = G⁻¹AᵀG, both pinned against a
   structurally-independent dense-transpose-plus-explicit-diagonal-G
   oracle (derivations/diagnostics/diag_p42_adjoint_oracle.py) — NOT a
   code-to-code comparison against another ORPHEUS adjoint path. This
   is the algebra-of-record ground for the equation.
.. vv-status: g-adjoint-definition documented

.. math::
   :label: g-adjoint-definition

   A^{\dagger} \;=\; G^{-1}\,A^{\mathsf T}\,G ,

acting on a composite ``bulk ⊕ boundary`` field (the timeless operator
carrier :class:`~orpheus.transport.full_field.FullField`; a driver's
history-bearing
:class:`~orpheus.transport.timed_full_field.TimedFullField` iterate reaches
``op.H`` via MRO, as in the reciprocity gate). Before
R5, ``op.H`` silently reduced to the plain **Euclidean** transpose
:math:`A^{\mathsf T}` because the SN operators advertised no
metric-bearing ``domain`` / ``codomain``; the wrapper had no metric to
read, so :math:`G` defaulted to the identity. R5 supplies the metric by
giving the FULL streaming leaf a direct-sum function space, and the
already-existing :class:`~orpheus.numerics.operator._AdjointOperator`
wrapper turns that into the correct G-adjoint **with no change to the
wrapper**.

.. admonition:: Key Facts (composite G-adjoint)
   :class: tip

   - **The Hilbert adjoint is the G-adjoint** :math:`A^{\dagger} =
     G^{-1} A^{\mathsf T} G`, NOT the Euclidean transpose. The
     reciprocity identity it satisfies is :math:`\langle A\psi,
     \varphi\rangle_G = \langle \psi, A^{\dagger}\varphi\rangle_G`.
   - **G is block-diagonal** on :math:`V = V_{\rm bulk}\oplus
     V_{\rm trace}`: bulk block :math:`G_{\rm bulk} = V_{\rm cell}\,w_n`
     (the full phase-space measure :math:`\mathrm dV\,\mathrm d\Omega`);
     trace block :math:`G_{\rm trace} = |\Omega\cdot\hat n_f|\,w_n` (the
     partial-current surface measure). Both carry :math:`w_n`; they
     differ only in the spatial measure (cell volume vs. oriented face).
   - **The carrier is** :class:`~orpheus.numerics.spaces.full_field_space.FullFieldSpace`
     — a direct sum that dispatches the metric **per block** to its
     leaf spaces (a pure composition, no new metric arithmetic).
   - **The metric-blind adjoint is unrepresentable — the metric lives
     on the shared space, NOT any leaf's domain.** Since P4.5 W-D, ``C``
     / ``S`` / ``F`` carry the SAME composite ``full_field_space`` as
     ``L`` / ``B`` (so the within-group
     :class:`~orpheus.numerics.operator.OperatorSum` guard *validates*
     the loss composition); the metric applies **once at the op level**
     because the :class:`~orpheus.numerics.operator._AdjointOperator`
     wrapper reads it off the *composite* ``domain`` / ``codomain`` of
     the summed operator, never per-leaf. Because every loss leaf now
     carries the composite metric, no composite adjoint is metric-blind,
     and a non-adjointable operand still makes the recursive
     :attr:`~orpheus.numerics.operator.LinearOperator.is_adjointable`
     report ``False``, so ``.H`` **raises**
     :class:`~orpheus.numerics.operator.MissingAdjoint` **eagerly at**
     ``.H`` **construction** — it never silently goes Euclidean. (The
     reachability table below predates the #112/#118/#276 adjoint work
     and is superseded; see its retraction note.)
   - **The trace metric is singular** on tangential ordinates
     (:math:`|\Omega\cdot\hat n| = 0`), so :math:`G^{-1}` is the
     **Moore–Penrose pseudo-inverse** (zero on the null space). This is
     **exact** for the adjoint: the inflow / outflow selectors exclude
     tangential ordinates, so those slots are identically zero in every
     matvec output.
   - **Discriminating gate:** the L11 wrong-metric control (drop
     :math:`|\Omega\cdot\hat n|` from the trace block) breaks
     reciprocity by :math:`6.4\times10^{-2}` (slab) /
     :math:`8.3\times10^{-1}` (sphere) / :math:`1.9\times10^{-3}` (cyl)
     — all :math:`\gg 10^{-3}`, proving the weighting is load-bearing.


The G-adjoint, derived from reciprocity
---------------------------------------

The Hilbert adjoint :math:`A^{\dagger}` of a linear operator on a space
with inner product :math:`\langle\cdot,\cdot\rangle_G` is *defined* by
the reciprocity (turn-over) identity

.. math::
   :label: g-adjoint-reciprocity

   \langle A\psi,\,\varphi\rangle_G
   \;=\;
   \langle \psi,\,A^{\dagger}\varphi\rangle_G
   \qquad\forall\,\psi,\varphi \in V .

The diagonal metric is represented by a (block-diagonal) Gram matrix
:math:`G`, so the inner product is :math:`\langle a,b\rangle_G = a^{\mathsf
T} G\, b`. Substituting into :eq:`g-adjoint-reciprocity` and solving for
:math:`A^{\dagger}` recovers :eq:`g-adjoint-definition` in three steps:

.. math::
   :label: g-adjoint-derivation

   \langle A\psi, \varphi\rangle_G
   = (A\psi)^{\mathsf T} G\, \varphi
   = \psi^{\mathsf T}\,(A^{\mathsf T} G)\,\varphi
   \;\overset{!}{=}\;
   \psi^{\mathsf T} G\, (A^{\dagger}\varphi)
   = \langle \psi, A^{\dagger}\varphi\rangle_G ,

so :math:`A^{\mathsf T} G = G\,A^{\dagger}` and therefore
:math:`A^{\dagger} = G^{-1} A^{\mathsf T} G`. When :math:`G = I` (the
Euclidean default) this collapses to the bare transpose
:math:`A^{\dagger} = A^{\mathsf T}` — which is precisely the
**wrong** adjoint for the SN composite, because the bulk and boundary
blocks are integrated against *non-uniform* measures (a cell of twice
the volume contributes twice the inner-product weight; a grazing
ordinate contributes :math:`|\Omega\cdot\hat n|` less surface current).
The whole point of R5 is to supply the correct :math:`G` so the
reciprocity holds under the *physical* measure, not the counting
measure of array indices.

This is the discrete twin of the continuous transport reciprocity
:math:`\int \mathrm dV\!\int \mathrm d\Omega\; \varphi\,(A\psi) = \int
\mathrm dV\!\int \mathrm d\Omega\; \psi\,(A^{\dagger}\varphi)` (Lewis &
Miller 1993, §3.7) — the phase-space integral :math:`\mathrm
dV\,\mathrm d\Omega` is what the bulk metric :math:`V_{\rm cell}\,w_n`
discretizes, and the surface partial-current integral :math:`\int_\Gamma
|\Omega\cdot\hat n|\,\mathrm dA\,\mathrm d\Omega` is what the trace
metric :math:`|\Omega\cdot\hat n_f|\,w_n` discretizes.


The block-diagonal metric — why each block carries its measure
--------------------------------------------------------------

The composite lives on the **direct sum** :math:`V = V_{\rm bulk}\oplus
V_{\rm trace}`, so its Gram matrix is **block-diagonal** (the bulk and
trace degrees of freedom are distinct coordinates — there is no
cross-term, the two integrals are over different domains):

.. math::
   :label: g-adjoint-block-metric

   G \;=\;
   \begin{pmatrix} G_{\rm bulk} & 0 \\ 0 & G_{\rm trace} \end{pmatrix},
   \qquad
   G_{\rm bulk} = V_{\rm cell}\,w_n ,
   \qquad
   G_{\rm trace} = |\Omega\cdot\hat n_f|\,w_n .

**The bulk block** :math:`G_{\rm bulk} = V_{\rm cell}\,w_n`. The bulk
inner product is the discretization of the full phase-space integral

.. math::

   \langle a, b\rangle_{G_{\rm bulk}}
   \;=\;
   \int_{\mathcal D}\!\mathrm dV
   \int_{4\pi}\!\mathrm d\Omega \; a\,b
   \;\longrightarrow\;
   \sum_{i\in\text{cells}} V_i
   \sum_{n\in\text{ordinates}} w_n \; a_{n,i}\,b_{n,i} .

The two quadratures factor: the **cell volume** :math:`V_i`
discretizes :math:`\mathrm dV` and the **angular quadrature weight**
:math:`w_n` discretizes :math:`\mathrm d\Omega`. The product
:math:`V_i\,w_n` is therefore the diagonal phase-space measure
:math:`\mathrm dV\,\mathrm d\Omega`. In code
(:meth:`SNMesh.full_field_space <orpheus.sn.mesh.augmented_mesh.SNMesh.full_field_space>`)
it is built as

.. code-block:: python

   g_bulk = w_n[:, None, None, None] * V[None, None, :, :]   # (N, 1, nx, ny)

— a ``(N, 1, nx, ny)`` array. The leading axis carries the per-ordinate
:math:`w_n`; the two spatial axes carry the per-cell :math:`V`; the
**energy-group axis is a singleton** because the phase-space measure is
**group-independent** (a group does not change a cell's volume or an
ordinate's solid angle). The singleton broadcasts over the energy axis
of the ``(N, ng, nx, ny)`` bulk tensor at metric-application time, so
the same :math:`(N,1,nx,ny)` weight serves every group with no
duplication. This is exactly the leading-axis broadcast convention of
:meth:`FunctionSpace._broadcast_metric <orpheus.numerics.space.FunctionSpace>`.

**The trace block** :math:`G_{\rm trace} = |\Omega\cdot\hat n_f|\,w_n`.
The boundary inner product is the discretization of the **partial-current
surface** integral

.. math::

   \langle a, b\rangle_{G_{\rm trace}}
   \;=\;
   \int_{\Gamma}\!\mathrm dA
   \int_{4\pi}\!\mathrm d\Omega \; |\Omega\cdot\hat n_f| \; a\,b
   \;\longrightarrow\;
   \sum_{f\in\text{faces}}
   \sum_{n} |\Omega_n\cdot\hat n_f|\,w_n \; a_{n,f}\,b_{n,f} .

The boundary degrees of freedom are angular fluxes *on a surface*, and
the physically meaningful surface functional is the **partial current**
:math:`J^{\pm} = \int_{\Omega\cdot\hat n \gtrless 0} |\Omega\cdot\hat
n|\,\psi\,\mathrm d\Omega`, so the surface measure carries the
**cosine factor** :math:`|\Omega\cdot\hat n_f|`. This is the same
:math:`|\Omega\cdot\hat n|`-weighted inner product under which the
reflective / white boundary operators are self-adjoint (see
:ref:`bc-extraction` and :mod:`orpheus.numerics.spaces.angular_trace_space`),
already populated in sub-step 4.1 (commit ``89b2f62``). R5 reuses it
verbatim as the trace block — it does **not** re-derive it.

**Both blocks carry** :math:`w_n` **— they differ only in the spatial
measure.** This is the structural symmetry that makes the direct-sum
metric a clean composition rather than two unrelated weightings: angular
integration is identical on both blocks (the same quadrature
discretizes :math:`\mathrm d\Omega` everywhere), and the spatial measure
specializes — a 3-D **volume** :math:`V_{\rm cell}` for bulk degrees of
freedom that live *in* a cell, a 2-D **oriented surface**
:math:`|\Omega\cdot\hat n_f|\,\mathrm dA` for trace degrees of freedom
that live *on* a face.


The carrier — :class:`FullFieldSpace`, a per-block metric dispatcher
--------------------------------------------------------------------

The metric :eq:`g-adjoint-block-metric` is carried by a new function
space :class:`~orpheus.numerics.spaces.full_field_space.FullFieldSpace`
(``orpheus/numerics/spaces/full_field_space.py``). It holds the two leaf
spaces — a bulk :class:`~orpheus.numerics.space.FunctionSpace` whose
``inner_product_weights`` is :math:`G_{\rm bulk}`, and the existing
:class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace` whose
``inner_product_weights`` is :math:`G_{\rm trace}` — and **overrides**
the three metric primitives to dispatch **per block**:

.. code-block:: python

   def apply_metric(self, x):            # G ⊙ x
       return self._rebuild(
           x,
           self.bulk_space.apply_metric(x.bulk.values),
           self.trace_space.apply_metric(x.boundary.values),
       )

Each override splits the composite field into its ``.bulk`` /
``.boundary`` blocks, routes each to its leaf space's metric method, and
rebuilds the composite. This is a **pure Pattern-2 composition**: the
direct-sum space owns only the *structure* (how to split and recombine);
it introduces **no new metric arithmetic** — the per-axis broadcast and
the pseudo-inverse already live in
:class:`~orpheus.numerics.space.FunctionSpace`. The inner product is
correspondingly the **sum of block inner products**
:math:`\langle a, b\rangle_G = \langle a_{\rm bulk}, b_{\rm
bulk}\rangle_{G_{\rm bulk}} + \langle a_{\rm trace}, b_{\rm
trace}\rangle_{G_{\rm trace}}` — the defining property of a direct-sum
Hilbert space.

The space is **duck-typed** on the composite field (``.bulk`` /
``.boundary`` leaves, each a frozen dataclass with a ``.values``
ndarray; rebuilt via :func:`dataclasses.replace`) so that the
``numerics`` layer never imports the ``transport`` layer — an
architectural firewall that keeps the operator-algebra primitives
domain-agnostic. The identity is the inherited ``(name, shape)`` tuple
(``name = "full_field"``, ``shape = (n_bulk + n_trace,)``), with the
block spaces as ``compare=False`` leaf metadata — so two composites over
meshes of the same total dimension compare equal and the
:class:`~orpheus.numerics.operator.OperatorSum` composition guard
accepts the full within-group loss ``L + C - S - B`` (every operand —
:math:`L`, :math:`C`, :math:`S`, and :math:`B` — reports the same
composite domain; P4.5 W-D gave the previously ``None``-spaced
:math:`C`/:math:`S`/:math:`F` real spaces and de-SN-ified the name from
``"sn_full_field"``).
The mesh exposes it as the cached property
:meth:`SNMesh.full_field_space <orpheus.sn.mesh.augmented_mesh.SNMesh.full_field_space>`.

**The wrapper is unchanged.** The whole apparatus plugs into the
**pre-existing** :class:`~orpheus.numerics.operator._AdjointOperator`,
which realizes ``A.H`` as

.. math::
   :label: g-adjoint-wrapper-action

   (A^{\dagger}\varphi)
   \;=\;
   \underbrace{G_V^{+}}_{\substack{\text{domain}\\\text{inverse-metric}}}
   \;\odot\;
   A^{\mathsf T}\!\Bigl(
   \underbrace{G_W}_{\substack{\text{codomain}\\\text{metric}}}\!\odot\,\varphi
   \Bigr) ,

calling ``codomain.apply_metric`` *before* the transpose and
``domain.apply_inverse_metric`` *after* (operator.py, ``_AdjointOperator.apply``).
The wrapper is metric-*representation*-agnostic: the SAME code path
serves the flat-ndarray spherical-harmonic :math:`(L+1, 2L+1)`
leading-axis metric AND the composite ``bulk ⊕ trace`` metric on a
structured field. R5 added a function space, not a line of adjoint
logic — the cleanest possible realization of the
:ref:`metric-lives-at-the-leaf <eigenvalue-posing>` principle.


The G-adjoint applies the metric once at the op level — the adjoint-axis predicate
----------------------------------------------------------------------------------

The subtle architectural point of R5 is **where** the metric lives in
the adjoint, and it survives the P4.5 W-D change to the bulk leaves'
domains. Since W-D, all five operators — :math:`L`
(:class:`~orpheus.sn.operators.streaming.StreamingOperator`), the collision
multiplier :math:`C = M[\sigma_t]`
(:class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`),
:class:`~orpheus.transport.operators.scattering.ScatteringOperator` (``S``),
:class:`~orpheus.transport.operators.fission.FissionOperator` (``F``), and
:math:`B` (:class:`~orpheus.sn.operators.boundary.SNBoundaryOperator`) —
carry the **same** composite ``full_field_space`` (threaded through
``from_solver_data`` / ``sn_mesh.full_field_space``), so the
within-group :class:`~orpheus.numerics.operator.OperatorSum` guard
*validates* the loss composition (``domain = None`` survives only on
the bare / test constructor). The architecturally interesting fact is
that this domain plumbing does **not** change where the metric is
applied in the adjoint, and that the metric-blind adjoint is made
**unrepresentable by the recursive** ``is_adjointable`` **predicate, not
by any leaf's domain**. Two mechanisms.

**(1) The metric applies ONCE at the op level — never per leaf.** The
:class:`~orpheus.numerics.operator._AdjointOperator` wrapper realizes
``op.H`` as :math:`G_V^{+}\odot(\cdot)^{\mathsf T}\!\bigl(G_W\odot(\cdot)\bigr)`,
reading :math:`G` off the **wrapped operator's** ``domain`` /
``codomain``. When the wrapped operator is the *sum* :math:`(L+C-B)`,
that is the sum's composite ``full_field_space`` (the
:class:`~orpheus.numerics.operator.OperatorSum` ``domain`` is the
first-non-``None`` of its operands — now redundant, since all operands
agree on the same composite). The adjoint therefore applies the metric
**once on the composite**:

.. math::
   :label: g-adjoint-sum-conjugation

   (L + C - B)^{\dagger}
   \;=\;
   G^{-1}\,(L + C - B)^{\mathsf T}\,G
   \;=\;
   G^{-1} L^{\mathsf T} G
   + G^{-1} C^{\mathsf T} G
   - G^{-1} B^{\mathsf T} G ,

distributing the **same** :math:`G^{-1}(\cdot)^{\mathsf T} G`
conjugation across every leaf in the sum. Although each leaf now
*advertises* the composite domain, the
:class:`~orpheus.numerics.operator._AdjointOperator` applies :math:`G`
at the **sum** level, never re-applying it per summand — the metric
weighting belongs to the *space*, applied once where the composite
enters and leaves, so a leaf carrying the composite domain is **not** a
double-application risk (the leaves' own ``apply_transpose``, where
defined, is the metric-blind Euclidean transpose; the metric is layered
on once by the sum's adjoint wrapper). The bulk leaves are pure
:math:`(N,ng,nx,ny)\to(N,ng,nx,ny)` endomorphisms whose transpose is
well-defined Euclidean-wise.

**(2) Every loss leaf carries the metric, so no composite adjoint is
metric-blind.** The concern the original design guarded against was an
adjoint taken over a sum that does **not** contain ``L`` and therefore,
under the *pre*-W-D design, carried no metric :math:`G`. That concern is
now **moot**: since P4.5 W-D every loss leaf — ``L``, ``C``, ``S``,
``F``, ``B`` — carries the SAME composite ``full_field_space``, whose
metric the :class:`~orpheus.numerics.operator._AdjointOperator` reads off
the composite domain. Every composite adjoint that constructs is
therefore metric-correct, whichever leaves it contains. The **surviving**
guard is general and lives on the adjoint axis: any operator with a
non-adjointable operand reports
:attr:`~orpheus.numerics.operator.LinearOperator.is_adjointable`
``= False``, and its :attr:`~orpheus.numerics.operator.LinearOperator.H`
raises :class:`~orpheus.numerics.operator.MissingAdjoint` **eagerly at
construction** — never silently Euclidean. This is
**illegal-states-unrepresentable** (Cardinal Rule 2) realized through the
recursive predicate, orthogonal to the W-D domain plumbing (which serves
the *forward* composition guard).

.. note:: **Supersession — the S/F adjointability update.** An earlier
   version of this section made the full prompt-loss adjoint
   :math:`(L + C - S - F - B)^{\dagger}` *intentionally unreachable* by
   having ``S`` and ``F`` carry **no** ``apply_transpose`` — a
   non-adjointable operand blocked the composite. **That mechanism is
   retired.** The #112 fission dyad-swap :math:`F^{\mathsf T}` and the
   #118 / #276 scattering Euclidean transpose :math:`S^{\mathsf T}`
   (via :attr:`~orpheus.transport.operators.scattering.ScatteringOperator.full_scatter_kernel`)
   gave both leaves a working ``apply_transpose``, so
   :class:`~orpheus.transport.operators.scattering.ScatteringOperator`
   and :class:`~orpheus.transport.operators.fission.FissionOperator` now
   report ``is_adjointable = True`` (``is_invertible`` still ``False`` —
   a source operator is not invertible). The metric-blindness concern is
   moot regardless, for two reasons: **(a)** every leaf carries the
   composite ``full_field_space`` metric (above); and **(b)** the
   within-group loss is never fused into a single
   :math:`(L + C - S - F - B)` operator —
   :func:`~orpheus.sn.coupled_system.build_within_group_system` returns the
   :class:`~orpheus.sn.coupled_system.WithinGroupSystem` record whose
   ``resolvent`` ``(L+C)`` and lagged ``gains`` ``(S, B_a)`` keep ``S`` /
   ``B`` as **lagged gains** and
   ``F`` handled at the eigenvalue / DSA **outer** layer (where the
   adjoint posing row daggers :math:`M = F` as :math:`M^{\dagger}`; see
   :ref:`eigenvalue-posing`), not via a within-group composite adjoint.


The singular trace metric and the pseudo-inverse — exactness
------------------------------------------------------------

The trace block :math:`G_{\rm trace} = |\Omega\cdot\hat n_f|\,w_n` is
**singular**: on a **tangential** ordinate (one with
:math:`\Omega\cdot\hat n_f = 0`, grazing the face) the cosine factor
vanishes, so that diagonal entry of :math:`G` is zero. A literal
:math:`G^{-1}` does not exist. R5 uses the **Moore–Penrose
pseudo-inverse** :math:`G^{+}` — :math:`1/G` where :math:`G \neq 0`, and
:math:`0` on the null space (:meth:`FunctionSpace.apply_inverse_metric
<orpheus.numerics.space.FunctionSpace>`):

.. code-block:: python

   nonzero = wb != 0.0
   return np.where(nonzero, x / np.where(nonzero, wb, 1.0), 0.0)

This is not an approximation — it is **exact** for the adjoint, by the
following argument. The pseudo-inverse zeroes the tangential components
of the adjoint output. That is the correct value because the tangential
trace slots are **identically zero in every matvec output** in the first
place: the boundary inflow / outflow selectors
(:meth:`AngularTraceSpace.outflow_indices_for_face
<orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace>`) classify an ordinate
as inflow (:math:`\Omega\cdot\hat n < -\epsilon`), outflow
(:math:`> +\epsilon`), or **tangential** (:math:`|\cdot|\le\epsilon`) —
and the boundary operators read/write **only** the inflow and outflow
slots. The tangential slots are never sourced. Consequently:

* the tangential components of :math:`A\psi` are zero (no operator
  writes them), so they carry **zero weight** in
  :math:`\langle\cdot,\cdot\rangle_G` anyway (:math:`G_{\rm trace} = 0`
  there);
* the pseudo-inverse returns zero on exactly those components, so the
  reconstructed adjoint output :math:`G^{+} A^{\mathsf T}(G\varphi)`
  agrees with the true G-adjoint on the **range** of the operator (the
  inflow ⊕ outflow subspace), and is zero on the orthogonal complement
  (the tangential null space) — which is where the true adjoint is also
  zero.

The pseudo-inverse and the true inverse-restricted-to-the-range
**coincide on the subspace the operator actually touches**, so the
reciprocity :eq:`g-adjoint-reciprocity` holds to round-off (it would
fail only if a matvec ever sourced a tangential slot — which the
selectors forbid). The trace-block residuals in the verification table
below (:math:`\le 3.6\times10^{-15}`) confirm there is no measurable
contamination from the null space.


Numerical evidence
------------------

The defining ground is the dense-probe oracle
:func:`validate_composite_adjoint <derivations.diagnostics.diag_p42_adjoint_oracle>`
(``derivations/diagnostics/diag_p42_adjoint_oracle.py``). It is
**structurally independent** of the production path: it assembles the
operator's dense matrix by probing :math:`op.\text{apply}` on unit
vectors, builds the diagonal metric :math:`G` **explicitly** from
:math:`V`, :math:`w_n`, and ``trace.inner_product_weights`` (it never
calls :meth:`FullFieldSpace.apply_metric
<orpheus.numerics.spaces.full_field_space.FullFieldSpace>`), and
compares ``op.H.apply(φ)`` against the explicit fold :math:`G^{+}\bigl(
op_{\rm dense}^{\mathsf T}\,(G\cdot\varphi)\bigr)` for
:math:`op = (L + C - B)`.

.. list-table:: ``op.H`` vs. the explicit block-diagonal G-fold :math:`G^{+}(op^{\mathsf T}(G\varphi))`, :math:`op = L + C - B`
   :header-rows: 1
   :widths: 20 22 22 36

   * - Geometry
     - bulk block :math:`|\Delta|_\infty`
     - trace block :math:`|\Delta|_\infty`
     - G-reciprocity :math:`\langle op\,\psi,\varphi\rangle_G = \langle\psi, op.H\,\varphi\rangle_G` (rel)
   * - slab (2-group)
     - :math:`7.1\times10^{-15}`
     - :math:`2.5\times10^{-16}`
     - :math:`6.5\times10^{-17}`
   * - sphere (2-group)
     - :math:`1.7\times10^{-13}`
     - :math:`3.6\times10^{-15}`
     - :math:`1.6\times10^{-15}`
   * - cylinder
     - :math:`2.8\times10^{-14}`
     - :math:`1.8\times10^{-15}`
     - :math:`6.8\times10^{-17}`

Both blocks match the explicit fold to round-off, and the defining
reciprocity holds to :math:`\le 1.6\times10^{-15}` across all three
geometries. Because the oracle's :math:`G` is built from raw mesh
quantities while production's :math:`G` is built by
:class:`FullFieldSpace`, the agreement also cross-validates the **metric
population** (a wrong :math:`V` or a transposed :math:`w_n` would shift
the bulk residual off round-off).

**The L11 wrong-metric negative control** is the discriminating gate
that proves the :math:`|\Omega\cdot\hat n|` weighting is **load-bearing**
and not an inert decoration. The control re-evaluates reciprocity under
a *deliberately wrong* trace metric — the angular weight :math:`w_n`
alone, with the :math:`|\Omega\cdot\hat n|` cosine factor dropped — while
``op.H`` is still the adjoint built for the **true** metric. A
correct-but-redundant weighting would leave reciprocity intact; a
load-bearing one must break it:

.. list-table:: L11 control — drop :math:`|\Omega\cdot\hat n|` from the trace metric (must break reciprocity, :math:`\gg 10^{-3}`)
   :header-rows: 1
   :widths: 30 35 35

   * - Geometry
     - reciprocity residual (rel)
     - verdict
   * - slab
     - :math:`6.4\times10^{-2}`
     - **broken** (:math:`\gg 10^{-3}`)
   * - sphere
     - :math:`8.3\times10^{-1}`
     - **broken** (:math:`\gg 10^{-3}`)
   * - cylinder
     - :math:`1.9\times10^{-3}`
     - **broken** (:math:`> 10^{-3}`)

All three break by orders of magnitude relative to the round-off
reciprocity of the correct metric — the cosine factor is doing real
work. (The cylinder margin is the smallest because its single curved
face has the narrowest spread of :math:`|\Omega\cdot\hat n|` over the
ordinate set; the slab and sphere carry the decisive controls, which is
why the L11 test gates on slab and sphere specifically.)

These results are pinned by two foundation-tagged test files:

* ``tests/sn/operators/test_g_adjoint_reciprocity.py`` — **12 tests**:
  the G-adjoint reciprocity on slab / sphere / cylinder / slab-2g /
  sphere-2g (5), a metric-population cross-check that
  ``op.codomain.inner_product`` matches an independent reference built
  directly from ``omega_dot_n`` / ``volumes`` (5), and the L11
  wrong-metric control on slab / sphere (2). The reciprocity inner
  products are evaluated with an **independent** Gram fold so a wrong
  *metric* cannot mask a wrong *adjoint*.
* ``tests/numerics/test_full_field_space.py`` — **6 tests** pinning the
  :class:`FullFieldSpace` identity semantics (flat direct-sum ``shape``,
  ``(name, shape)``-only identity with ``compare=False`` block
  metadata, dict-key usability, the composite's own
  ``inner_product_weights`` staying ``None``).

Both files are ``@pytest.mark.foundation`` (software invariants over the
operator-algebra ground truth — they carry **no** ``verifies()`` label
because there is no solver-level theory equation being checked; the
defining identity is the algebra of :eq:`g-adjoint-definition` itself,
anchored to the structurally-independent oracle, not to a discretization
claim). The **forward path is bit-identical**: the diamond-difference
regression suite (69 passed) confirms the adjoint addition does not
perturb the forward matvec — R5 added a *new* capability (``.H``), it
did not touch ``apply``.

.. note::

   This section discharges the **adjoint-metric** half of the "what
   remains for O.2" list in :ref:`bc-extraction-operator-output-o2`
   together with **Gate-1.3** (the O.2 adjoint verification gate). The
   remaining open item is the **residual column** —
   :meth:`AngularBoundaryResidual.from_balance
   <orpheus.transport.residuals.angular_boundary_residual.AngularBoundaryResidual.from_balance>`
   has no operator-output consumer until the O.2 named-composition
   driver types the affine boundary balance at the solver level. The
   G-adjoint of this section is what gives the
   :ref:`adjoint posing row <eigenvalue-posing>` its
   :math:`G`-weighted transpose for free — the daggered eigenvalue
   problem :math:`A_{\rm loss}^{\dagger}\psi^{\dagger} = \lambda
   M^{\dagger}\psi^{\dagger}` now has a concrete, verified ``.H`` to
   build on.


.. _wavefront-flux-cochain:

The interior face-flux cochain — :math:`C^1_{\rm int}`
======================================================

.. note:: **Succession (S6.4(f), Issue #222, 2026-06-10) — the typed
   carrier retired; the concept survives in its two native
   realizations.**

   The interior face cochain :math:`C^1_{\rm int}` was, from #205 Phase 5
   through S6.4, carried by a dedicated named field ``WavefrontFlux``
   (on the space ``InteriorFaceSpace``). **That type is retired** — the
   modules ``orpheus/transport/fields/wavefront_flux.py`` and
   ``orpheus/numerics/spaces/interior_face_space.py``, and the 25
   foundation tests in ``tests/transport/fields/test_wavefront_flux.py``,
   are deleted. The cochain **mathematics** below — the biproduct
   :math:`C^1 = C^1_{\rm int} \oplus C^1_\partial`, the trace algebra
   :math:`\iota_*` / :math:`\iota^*`, the flux-only single-role
   rationale, the storage × role × locus grid — **remains valid theory**.
   Only the Python carrier is gone, and the historical derivation is
   kept below in past tense because it explains *why* the cochain frame
   is the right one.

   **Why it retired.** The S6.4 walk re-layering (see
   :ref:`sweep-dispatch-relayering` in :doc:`discrete_ordinates`)
   dissolved the type's two load-bearing verbs — the whole-trace
   :math:`\iota_*` ``seed`` (read the octant inflow into the
   domain-edge slots) and the whole-trace :math:`\iota^*` ``absorb``
   (capture the domain-edge outflow) — into the shared ``_OctantWalk``
   frame, where the per-octant inflow read, the per-octant outflow
   shed, and the **single** :ref:`O.4b <bc-extraction>` boundary block
   live once for every walk. The type's whole-:math:`N` **mesh-bound**
   storage then had no place in the **per-octant** walk: an octant
   transient cannot be a persistent mesh-bound field.

   **Where the concept lives now** — two native realizations, each
   truer to a different facet of what the type held:

   * **The values AT the moving wavefront** —
     ``orpheus.sn.loss_representation.sweep_graph._MovingFrontier``, the rolling
     :math:`(d{-}1)`-frontier (per-level ``seed`` /
     ``shed``). Arguably the *truer* realization: the retired type
     held the wavefront's complete **history**; the frontier holds the
     **front itself**, ping-ponged by parity across two active levels.
   * **The full cochain history (the oracle's fuller view)** —
     ``FullFieldWavefront._octant_face_cochain`` raw per-axis buffers
     in ``orpheus.sn.loss_representation`` (the in-edge :math:`\iota_*`
     seed; ``_edge_outflow`` is the :math:`\iota^*` extraction),
     retained for verification cross-checks of the production
     ``_MovingFrontier`` window.

   The type was *a useful concept while it lived* — it named the
   :math:`\iota_*` / :math:`\iota^*` trace operators the raw-numpy
   sweep applied by hand, and it is the substrate on which the
   biproduct decomposition was first articulated. The succession keeps
   that articulation; only the now-redundant carrier is dropped.

Wave O step #205 Phase 5 (`Issue #208
<https://github.com/deOliveira-R/ORPHEUS/issues/208>`_, commits
``478723d`` mint / ``992b0c0`` 2-D sweep / ``0e3e16c`` 2-D matvec,
2026-06-04) typed the SN wavefront sweep's **interior** cell-face
angular fluxes — historically the raw ephemeral numpy arrays ``psi_x``
``(N, ng, nx{+}1, ny)`` / ``psi_y`` ``(N, ng, nx, ny{+}1)`` — as a
named field ``WavefrontFlux``.
This was the **face-locus** sibling of the boundary-block typing of
:ref:`bc-extraction` (cell + trace) and the operator-output typing of
:ref:`bc-extraction-operator-output-typing`: where those typed the
*cell* state and the *operator outputs*, this typed the *interior
face* state that the wavefront sweep propagated between them. It killed
the ``coding-elegance`` Pattern-3 anti-pattern (a flux-bearing tensor
with no type identity) and **named the trace operator the sweep
applied by hand** — the algebra that the S6.4 walk frame later
absorbed.


The native frame — discrete exterior calculus / cochains
--------------------------------------------------------

The native mathematical structure of the SN face fluxes (validated by
the cross-domain-attacker, frame memo
``field_role_typing_faceflux_frames.md``) is **discrete exterior
calculus**. The per-ordinate angular flux crossing faces is a primal
**1-cochain** — an assignment of a value to each oriented face:

.. (vv-status rationale) Structural framing of the SN face fluxes as a
   primal 1-cochain. This is a definitional / representational identity
   (the named-field typing — now carried by _MovingFrontier /
   _octant_face_cochain after the WavefrontFlux carrier retired at
   S6.4(f)), not a solver claim; the verifiable content is the
   bit-identity of the typed walk against the raw psi_x/psi_y walk,
   pinned by the octant-equivalence and Gate-K suites below (the 25
   foundation tests of the retired test_wavefront_flux.py went with the
   type; the window ≡ full oracle now pins the cochain walk).
.. vv-status: wavefront-cochain-primal documented

.. math::
   :label: wavefront-cochain-primal

   \psi^{(1)}_\Omega \;\in\; C^1
   \qquad (\text{a value per oriented face, per ordinate } \Omega).

The cell-average :math:`\bar\psi` is a **0-cochain**
:math:`\bar\psi^{(0)} \in C^0`, and the diamond-difference closure
:math:`\psi^{\rm out} = 2\bar\psi - \psi^{\rm in}` is the discrete
statement that :math:`\bar\psi` is the *midpoint* of the face-pair
bounding the cell — i.e. the **averaging map**
:math:`C^1 \to C^0` co-restricted to the cell (the lowest-order
Whitney 1-form → 0-form reduction, the trapezoidal/midpoint Hodge star
on a 1-D edge). The :ref:`balance-cartesian-2d` DD already realises
this; the cochain frame simply names the spaces.


The biproduct :math:`C^1 = C^1_{\rm int} \oplus C^1_\partial`
-------------------------------------------------------------

The full face cochain splits as a **biproduct** into the interior and
the boundary 1-chains:

.. (vv-status rationale) The face-cochain biproduct decomposition. A
   structural identity mirroring the cell+trace direct sum of
   :eq:`bc-extraction-direct-sum-state` one locus down; the verifiable
   content is the round-trip + projection laws below, pinned by
   test_absorption_is_identity / test_pi_int_after_injection_is_zero_2d.
.. vv-status: wavefront-cochain-biproduct documented

.. math::
   :label: wavefront-cochain-biproduct

   C^1 \;=\; C^1_{\rm int} \;\oplus\; C^1_\partial,
   \qquad
   \begin{cases}
     C^1_{\rm int} & (\text{interior faces}), \\
     C^1_\partial  & (\text{domain-edge faces, } \texttt{AngularBoundaryFlux}),
   \end{cases}

coupled by the injection / projection at the domain edges. This is the
**same biproduct** as the cell+trace direct-sum transport state
:eq:`bc-extraction-direct-sum-state`
(:math:`V = V_{\rm bulk} \oplus V_{\rm boundary}`), realised **one
locus down** at the face level. The two loci nest:

.. list-table:: The two biproduct loci (cell+trace ‖ face)
   :header-rows: 1
   :widths: 22 30 26 22

   * - Locus
     - Interior summand
     - Boundary summand
     - Coupling
   * - **cell + trace**
       (:eq:`bc-extraction-direct-sum-state`)
     - :math:`V_{\rm bulk}`
       (:class:`~orpheus.transport.fields.angular_flux.AngularFlux`)
     - :math:`V_{\rm inflow} \oplus V_{\rm outflow}`
       (:class:`~orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux`)
     - the streaming sweep + sibling :math:`-B`
   * - **face**
       (:eq:`wavefront-cochain-biproduct`)
     - :math:`C^1_{\rm int}`
       (``_MovingFrontier`` front / ``_octant_face_cochain`` history)
     - :math:`C^1_\partial`
       (:class:`~orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux`)
     - the trace operators :math:`\iota_*` / :math:`\iota^*`

:class:`AngularBoundaryFlux` is the boundary summand at **both** loci — it is
literally the domain-edge faces, which are exactly the
boundary trace of the cell+trace state. The interior summand differs:
the cell biproduct carries the cell-centre flux, the face biproduct
carries the interior *face* flux. The boundary persists (it carries
the SI / Krylov iterate across calls); the interior is **ephemeral**
(rebuilt each sweep — which is *why*, at S6.4(f), no mesh-bound type
survived for it), so the two summands have different lifetimes —
:eq:`wavefront-cochain-biproduct` is therefore a lifetime split as
well as a spatial split.


The trace operators :math:`\iota_*` (seed) / :math:`\iota^*` (absorb)
---------------------------------------------------------------------

The boundary trace is the **pullback** :math:`\iota^*` of the interior
cochain under the inclusion :math:`\iota \colon \partial\Omega
\hookrightarrow \Omega`. In the discrete setting :math:`\iota^*` is
"select the domain-edge faces"; its adjoint injection :math:`\iota_*`
is "write the boundary trace into the domain-edge faces". These are
exactly what the pre-typing sweep did by hand —
``psi_x[:, :, 0, :] = boundary.face_view("xmin")`` (the seed) and the
write-back (the absorb). The retired ``WavefrontFlux`` named them as
``seed`` / ``absorb`` methods; after S6.4(f) the same verbs live on the
two native realizations — ``_MovingFrontier.seed`` / ``.shed`` (the
per-level front) and ``FullFieldWavefront._octant_face_cochain`` (the
in-edge :math:`\iota_*` seed) / ``_edge_outflow`` (the :math:`\iota^*`
extraction):

.. list-table:: The trace operators (now on the realizations)
   :header-rows: 1
   :widths: 14 22 28 36

   * - Symbol
     - Realization verb
     - Direction
     - Role in the biproduct
   * - :math:`\iota_*`
     - ``_MovingFrontier.seed`` /
       ``FullFieldWavefront._octant_face_cochain`` in-edge
     - :math:`C^1_\partial.\text{inflow} \to C^1_{\rm int}`
       domain-edge slots
     - the **injection** :math:`\iota_\partial` of the biproduct
   * - :math:`\iota^*`
     - ``_MovingFrontier.shed`` /
       ``FullFieldWavefront._edge_outflow``
     - :math:`C^1_{\rm int}` domain-edge slots
       :math:`\to C^1_\partial.\text{outflow}`
     - the **projection** :math:`\pi_\partial` of the biproduct

Historically the retired type routed both through a single
single-source-of-truth ``_edge_slot`` face-to-edge map, so the
injection and projection could not desync, and exposed a third read
accessor ``edge_view`` that returned the domain-edge slot as a
zero-copy view *without* copying into a :class:`AngularBoundaryFlux` — the
2-D matvec used it to difference ``edge_view(face) − given`` when
emitting its active-trace boundary residual (so there was no hardcoded
``psi_x[:, :, 0, :]`` literal at the call site). After S6.4(f) the same
single-map discipline is the shared ``_OctantWalk`` frame's: the
per-octant inflow read and outflow shed are the ONE boundary block
every walk uses, so the injection/projection cannot desync there
either.

The two biproduct laws follow, and are **provable** rather than
coincidental:

.. (vv-status rationale) The two biproduct identities — absorption =
   identity (project-after-inject) and projection-annihilates-the-
   strictly-interior (the off-diagonal block is zero). Structural laws
   of the biproduct; were pinned by test_absorption_is_identity (slab /
   sphere / 2-D box) and test_pi_int_after_injection_is_zero_2d in the
   retired test_wavefront_flux.py. After S6.4(f) the same two laws are
   the content of the window ≡ full-field oracle (the seed/shed
   round-trip equals the full-cochain seed/extract, bit-identically).
.. vv-status: wavefront-cochain-biproduct-laws documented

.. math::
   :label: wavefront-cochain-biproduct-laws

   \iota^* \circ \iota_* \;=\; \mathrm{id}
   \quad (\text{on the boundary chain}),
   \qquad
   \pi_{\rm int} \circ \iota_\partial \;=\; 0.

The first — :math:`\iota^* \circ \iota_* = \mathrm{id}` — IS the
"absorption = identity" fact (``seed`` then ``absorb``, with no
wavefront walk between, returns the boundary trace unchanged). It was
an *observation* under the raw-numpy seed/absorb, then a **named
biproduct law** under the typed ``WavefrontFlux``; after S6.4(f) it is
the seed/shed round-trip the ``_OctantWalk`` frame performs once per
octant. The second — :math:`\pi_{\rm int} \circ
\iota_\partial = 0` — is the biproduct's off-diagonal-zero condition:
injecting the boundary leaves the **strictly**-interior faces
(positions :math:`1 \ldots n{-}1` along each axis) untouched at zero
(in the full-cochain realization, ``_octant_face_cochain`` zero-inits
every interior slot and seeds only the in-edge, so the strictly-interior
faces are untouched by construction).


Why the interior cochain is flux-only (single role)
---------------------------------------------------

The interior cochain carries the **flux** role only — unlike the
boundary trace, it has no source / residual leaves (the retired
``WavefrontFlux`` was, accordingly, flux-only; the surviving
``_MovingFrontier`` / ``_octant_face_cochain`` buffers are likewise
plain flux arrays). The reason is
structural, and it explains the role grid of
:ref:`bc-extraction-operator-output-typing` from the cochain side. Per
the second native frame (sparse triangular factorisation), the interior
face flux is the **off-diagonal coupling of the per-octant
lower-triangular streaming factor** :math:`L_{\rm oct}`: the entry
coupling a cell to its upwind neighbour is the flux on the shared
interior face. That factor is **re-formed each sweep** (forward
substitution down the upwind DAG), so the interior face flux is a
*transient of the factorisation*, not a persistent state with a defect
of its own. The role grid (flux / source / residual) is a **0-cochain
(cell) concept**: a residual arises from a balance of cell-level
reaction rates against a cell-level source. Only the **boundary
1-chain** inherits the role grid — and only because the boundary
consistency residual

.. math::

   r_\Gamma \;=\; \psi.\text{inflow}
                  \;-\; B\,\psi.\text{outflow}
                  \;-\; q.\text{inflow}

(:eq:`bc-extraction-two-residuals`) lives on the boundary faces. The
strictly-interior faces carry no such balance, so the interior cochain
is flux-only by construction — ``illegal-states-unrepresentable``: there
was no ``InteriorFaceResidual`` to mistype an interior face as (and
under the surviving plain-array realizations there is no role grid to
violate at all).


The storage × role × locus grid (Issue #205), extended with the face locus
---------------------------------------------------------------------------

`Issue #205
<https://github.com/deOliveira-R/ORPHEUS/issues/205>`_ frames the
cross-method field architecture as a **storage × role × locus**
classification: a typed transport field is pinned by *where* it lives
(the **locus** — cell, face, boundary), *what role* it plays (flux /
source / residual), and *how it is stored* (persistent vs ephemeral).
Phase 5 adds the **face** locus row:

.. list-table:: Storage × role × locus — the typed field grid (#205)
   :header-rows: 1
   :widths: 16 18 22 22 22

   * - Locus
     - Storage
     - Flux role
     - Source/sink role
     - Residual role
   * - **cell**
       (0-cochain :math:`C^0`)
     - persistent (the iterate)
     - :class:`~orpheus.transport.fields.angular_flux.AngularFlux`
       / :class:`~orpheus.transport.fields.scalar_flux.ScalarFlux`
     - :class:`~orpheus.transport.source_sinks.angular_source_sink.AngularSourceSink`
     - :class:`~orpheus.transport.residuals.angular_residual.AngularResidual`
       (minted, O.2 consumer)
   * - **face — interior**
       (1-cochain :math:`C^1_{\rm int}`)
     - **ephemeral** (rebuilt each sweep)
     - ``_MovingFrontier`` front /
       ``_octant_face_cochain`` history
       (**#205 Phase 5** ``WavefrontFlux``, retired S6.4(f))
     - — (off-diagonal of :math:`L_{\rm oct}`, no source role)
     - — (no cell-balance defect on interior faces)
   * - **face — boundary**
       (1-cochain :math:`C^1_\partial`)
     - persistent (carries the iterate trace)
     - :class:`~orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux`
     - :class:`~orpheus.transport.source_sinks.angular_boundary_source_sink.AngularBoundarySourceSink`
       (B.5.2)
     - :class:`~orpheus.transport.residuals.angular_boundary_residual.AngularBoundaryResidual`
       (minted, O.2 consumer)

The grid makes the flux-only-single-role rationale (above) visible: the
interior-face row has only the flux cell populated, because the
off-diagonal-of-:math:`L_{\rm oct}` reading forbids a source role and
the no-cell-balance reading forbids a residual role. The boundary-face
row is fully populated only because the boundary 1-chain hosts the BC
consistency residual. The face locus mirrors the cell locus on storage
(interior-face ephemeral ‖ boundary-face persistent is the face-level
analogue of the bulk iterate persisting while operator outputs are
transient).


Field + views, NOT per-face objects
------------------------------------

The retired ``WavefrontFlux`` stored a **single flat backing buffer**
(``space.layout.total_size``); the per-axis face fields were zero-copy
reshape views (a ``face(axis)`` accessor).
The cross-domain-attacker **rejected** a per-face Python object on three
independent grounds — a rejection that **still binds** the surviving
realizations (``_MovingFrontier`` and ``_octant_face_cochain`` are both
plain dense arrays indexed by the batch, never per-face objects):

1. **Vectorisation (load-bearing).** The unit of operation is the
   ``(N_oct, ng, n_diag)`` wavefront batch; a per-face object would
   reintroduce the per-cell / per-face Python fold that caused a
   10–20× slab regression in earlier work. The hot-path walk indexes
   the ``face(axis)`` views with byte-identical fancy-indexing to the
   legacy raw ``psi_x`` / ``psi_y``.
2. **The cochain frame is storage-granularity-indifferent.** A cochain
   is an assignment of values to faces; whether stored as one flat
   buffer or per-face objects is an implementation choice the frame
   does not constrain — so the vectorisation argument settles it for
   the flat buffer.
3. **Biproduct consistency.** The
   :eq:`wavefront-cochain-biproduct-laws` projections are clean
   slice-views only if both summands live on the same flat-buffer
   substrate.

This mirrored :class:`AngularBoundaryFlux`, which still uses
``flat-buffer + FaceLayout + face_view``. The interior space was the
retired ``InteriorFaceSpace`` — the **layout-on-space** (A.5) sibling
of :class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace` minus the
``omega_dot_n`` directional selectors (the interior cochain has no
inflow/outflow partition — it is flux-only). It was axis-parametric: its
``interior_layout`` builder took the axis count as a parameter (one
face-normal field per active axis). The dimension-genericity that this
parametricity foreshadowed is now carried *directly* by the
realizations: ``_MovingFrontier`` is built from a mesh-time
``_FrontierPlan`` whose ``(d{-}1)``-frontier slab is constructed for any
``d`` (a point at d=1, a line at d=2, a surface at d=3), and
``FullFieldWavefront._octant_face_cochain`` allocates one ``n_a + 1``
buffer per axis for any ``ndim`` — so a future 3-D Cartesian wavefront
sweep is a *parameter*, not a new code path
(``feedback_unify_after_two_instances``: 1-D + 2-D are two working
instances, 3-D is the validating third).


.. _wavefront-flux-honest-scope:

Honest scope — a representation win, NOT a speed/rate/parallelism win
---------------------------------------------------------------------

.. warning::

   The interior-cochain **typing** was a **representation / elegance win
   only**. It did **NOT**:

   * change the asymptotic cost of the sweep,
   * recover the source-iteration convergence rate,
   * enable parallelism on its own.

   The cross-domain-attacker stated this as an *absence*, not a hedge:
   the wavefront DAG already gives the optimal sweep schedule, and the
   typing relocated no arithmetic. The seed/absorb stays an **inherent
   cheap** :math:`O(\text{boundary faces})` copy — negligible against
   the :math:`O(\text{volume})` sweep — at the **persistent-boundary /
   ephemeral-interior lifetime split**. True zero-copy between the two
   summands is *precluded* (and unnecessary) because
   :class:`AngularBoundaryFlux` persists across SI iterations while the
   interior cochain is rebuilt each sweep — the very lifetime split that
   later left it with **no** mesh-bound type after S6.4(f).

   The asymptotic-cost / peak-memory wins are sought elsewhere: the
   **persistent** SI iterate is shrunk by the orthogonal
   :ref:`angular windowing <sn-angular-windowing>` of Phase 5a (the
   iterate lives in moment space, :math:`N \to (L{+}1)(2L{+}1)`), and the
   **per-sweep transient** interior cochain
   itself — the dominant peak-memory cost noted here — was the target of
   Phase 5b's storage-B rolling moving-frontier window, which never
   materializes the whole interior cochain at once and is the 3-D
   enabler. That window IS the ``_MovingFrontier`` that the interior
   cochain now lives in: the production realization *is* the
   memory-win, with the full-cochain ``_octant_face_cochain`` kept only
   as the verification oracle.

The payoff *while the type lived* was the **type**: a named field,
typed :math:`\iota_*` / :math:`\iota^*`, code that reads like the
cochain math, and illegal-states-unrepresentable (the flux-only
constraint of the interior locus enforced by there being no
interior-face residual leaf). It was also framed as the **clean
substrate** the SI Gauss–Seidel recovery would land on: the
``(octant × face)`` reflective-graph G-S schedule as an explicit
composition (``sweep octant → ι* absorb → −B reflect → ι_* seed next
octant``) rather than an implicit buffer-timing dance. The S6.4 walk
re-layering delivered exactly that explicit composition — but as the
shared ``_OctantWalk`` frame's per-octant ``read inflow → shed outflow
→ −B boundary block`` sequence, **without** a standalone interior-face
type. The composition survived; the carrier did not. That recovery
remains where the actual convergence-rate win is sought (a separate,
research-tagged step).


Numerical evidence — type-only ⟹ bit-identical
-----------------------------------------------

Like the B.5.2 retype, Phase 5 wraps the **same** buffers: the
``face(axis)`` views share memory with the flat backing, and the
wavefront walk indexes them with byte-identical fancy-indexing, so
``.values`` are unchanged — only the type and the named
:math:`\iota_*` / :math:`\iota^*` differ. The change is verified by
**bit-identity (``np.array_equal``) by inheritance** from the
already-verified raw-numpy path:

.. list-table:: Phase 5 bit-identity gates
   :header-rows: 1
   :widths: 30 38 32

   * - Gate
     - Pins
     - Evidence
   * - **Phase 0 de-risk**
       (``diag_phase0_wavefront_derisk.py``, diagnostic)
     - typed :math:`\iota_*` / :math:`\iota^*` round-trip + full bare
       sweep on a flat-buffer backing
     - :math:`\max|\Delta| = 0.0` on angular / scalar / boundary,
       all reflective + vacuum cases; ``shares_memory`` zero-copy;
       transposed-seed negative control breaks (non-vacuous)
   * - **2-D octant equivalence**
       (``test_2d_octant_sweep_equivalence.py``)
     - the typed 2-D sweep reproduces the legacy octant snapshots
     - bit-identical octant snapshots
   * - **Gate-K** :math:`k_\infty = 1.875`
       (``test_2d_l2_matvec_correctness``)
     - the 2-D matvec recovers the closed-form homogeneous eigenvalue
       (:math:`k_\infty = \nu\Sigma_f / \Sigma_a`, closed-form pillar)
     - bit-identical to the pre-typing value
   * - **matvec ≡ sweep** + **BC extraction**
       (``test_bc_extraction_2d`` / ``_matvec``)
     - the typed matvec and sweep are ONE discretisation; the bare-2-D
       BC block matches
     - 126 passed
   * - **foundation suite**
       (``tests/transport/fields/test_wavefront_flux.py`` — RETIRED
       at S6.4(f) with the type)
     - units / class identity / field+views / the two biproduct laws /
       the round-trip pin + L11 negative control / axis-parametricity
       (1-D / 2-D / 3-D one path)
     - was all ``@pytest.mark.foundation`` (software invariants, no
       theory ``:label:``); the biproduct-law content these pinned now
       lives in the ``window ≡ full-field`` oracle
   * - **L16 perf**
       (``diag_l16_wavefront_microbench.py``, diagnostic)
     - no per-cell / per-face Python crept into the hot path
     - typed / raw median wall-clock ratio :math:`\approx 1.00`
       (within the +5 % gate)

The closed-form Gate-K eigenvalue is the only *correctness* (not
merely *equivalence*) anchor here: it is verified against
:math:`k_\infty = \nu\Sigma_f / \Sigma_a` (the closed-form pillar — MMS
does **not** prove eigenvalues). The remaining gates are bit-identity
by inheritance (``vv-principles`` § "Bit-identity vs
principled-equivalence", criterion: the implementation is unchanged —
the wrapper wraps the same buffers, so the values are byte-identical,
not merely close). The **A2D-1 source-hash** deliberate-edit tripwire
was refreshed (``f683f229…`` → ``12697ab3…``, behaviour-neutral — the
hash pins that any future edit to the 2-D matvec source is intentional).


The 1-D sweep is a scan, not a wavefront — deferred to ``nd_foundation``
------------------------------------------------------------------------

.. note::

   The retired ``WavefrontFlux`` typed the **2-D** wavefront sweep +
   matvec only. The **1-D** sweep is a parallel-prefix **scan** (a
   different fold): its interior fluxes are transient chain-ordered scan
   output (``(nx, K, ng)``, no persistent interior buffer at all), the
   structural antithesis of the cochain ``(N, ng, nx{+}1)``. Forcing the
   type into the 1-D scan would have been a wrong-fit (it would risk the
   L16 ``cumprod`` efficiency and multiply concepts), so the 1-D
   *implementation* kept its scan fold. The one shared seam (the
   boundary-trace exchange + DD-closure averaging) was deferred to a
   future ``nd_foundation`` session that would re-express **both** folds
   as one :math:`d`-generic walk, under the hard constraint that the
   collapse must not regress the :math:`d=1` scan's parallel-prefix
   efficiency. **That collapse has since landed** (the ``CumprodScan``
   d=1 sibling and the ``_DAGWavefront`` family share the per-octant
   DAG; see :ref:`sweep-dispatch-relayering` in
   :doc:`discrete_ordinates`), which is part of *why* the standalone
   ``WavefrontFlux`` carrier became redundant at S6.4(f): the
   :math:`d`-generic walk supplies the front (``_MovingFrontier``) and
   the oracle history (``_octant_face_cochain``) for every :math:`d`. A
   future 3-D *Cartesian* sweep IS a wavefront and rides the same
   ``_MovingFrontier`` window directly (the :math:`(d{-}1)`-frontier
   slab is built for any :math:`d`).


.. _sn-angular-windowing:

Angular windowing — the SI iterate lives in moment space (Phase 5a)
===================================================================

Wave O step #205 **Phase 5a** (commits ``93807aa`` factoring / ``b97d4f9``
eigenvalue inner / ``13ca001`` fixed-source inner, 2026-06-07) is a
**moment-reduction** of the SN within-group source-iteration *iterate*.
It is **orthogonal** to the :ref:`interior face-flux cochain <wavefront-flux-cochain>`
above: where that types the *per-ordinate* interior face flux a single
sweep propagates (and explicitly frames the interior cochain as
per-ordinate, :math:`\psi^{(1)}_\Omega \in C^1`), Phase 5a observes that
the **persistent** iterate the source iteration carries *between* sweeps
does not need all :math:`N` ordinates — it needs only the
spherical-harmonic moments the scattering operator consumes. The
held iterate's angular dimension drops :math:`N \to (L+1)(2L+1)`, and
the iterate becomes :class:`~orpheus.transport.fields.harmonic_moment_flux.HarmonicMomentFlux`
instead of :class:`~orpheus.transport.fields.angular_flux.AngularFlux`.

.. admonition:: Key Facts (angular windowing)
   :class: tip

   - **The SI fixed point lives in moment space.** Within-group source
     iteration is :math:`\psi_{k+1} = (L{+}C)^{-1}(S\,\psi_k + B\,\psi_k
     + q)`. The scattering source :math:`S\,\psi` is a pure function of
     the flux moments :math:`\phi_\ell^m = (M\psi)_\ell^m` — the
     per-ordinate iterate :math:`\psi` carries strictly more than the
     iteration consumes.
   - **Hold the moments, not the ordinates.** The persistent iterate is
     the moment tensor :math:`\phi \in (L{+}1, 2L{+}1, n_g, n_x, n_y)`,
     not :math:`\psi \in (N, n_g, n_x, n_y)`. Measured **18.3×**
     persistent-iterate shrink at :math:`N = 110`, :math:`L = 1`
     (:math:`N / (L{+}1)(2L{+}1) = 110/6`).
   - **The per-step source is bit-identical.** :math:`S` consuming the
     moments equals :math:`S` consuming the full angular flux **bit for
     bit** (0 ULP) under the ORPHEUS unnormalized-harmonic convention.
     The only non-bit-identical change is the SI *convergence test*,
     which moves to the moment :math:`L^2` — *more* principled, not a
     regression.
   - **2-D Cartesian only (load-bearing).** Windowing is valid where the
     sweep is a **direct** solve with no per-ordinate-iterate seed.
     Curvilinear (1-D sphere / cylinder) **must** stay full-angular: the
     Morel–Montry Carlson coupled-pole closure seeds from the previous
     iterate's per-ordinate :math:`\psi` at :math:`\mu = -1`, which the
     moment tensor cannot reconstruct. The Krylov path stays
     full-angular too. Gated on the genuine ``sn_mesh.is_cartesian and
     sn_mesh.ndim == 2`` (C5.4 / #225 — the earlier ``reduced is None``
     proxy was also true at 3-D Cartesian).
   - **Interior-bulk only.** The reflective :math:`B` coupling reads the
     full per-ordinate boundary *trace*; windowing reduces only the
     interior bulk. The biproduct :eq:`wavefront-cochain-biproduct`
     keeps the trace :math:`C^1_\partial` a distinct, **un-reduced**
     summand.
   - **A representation + typed-state win, NOT yet a peak-memory win.**
     5a shrinks the *persistent* iterate (18.3×) and makes its type
     honest. The *peak* memory drops only modestly (~1.2× measured)
     because the per-sweep **transient** full-angular machinery still
     dominates — that transient is the target of Phase 5b (interior-face
     cochain) and Phase 5c
     (:ref:`full-angular output <sn-angular-windowing-in-sweep-accumulation>`,
     the 3.06× linear peak win).


.. _sn-angular-windowing-fixed-point:

The within-group fixed point lives in moment space
--------------------------------------------------

The within-group source iteration solves, for each outer step, the
fixed-point problem

.. (vv-status rationale) The within-group source-iteration fixed point.
   This is the governing iteration the windowing reorganizes; the
   verifiable content is the SI ≡ Krylov cross-check (Krylov stays
   full-angular) and the closed-form k_inf eigenvalue, not the rendered
   recurrence. Documented, not orphan-gated.
.. vv-status: si-within-group-fixed-point documented

.. math::
   :label: si-within-group-fixed-point

   \psi_{k+1}
   \;=\; (L + C)^{-1}\!\left( S\,\psi_k + B\,\psi_k + q \right),

where :math:`L + C` is the within-group **invertible resolvent** (the
streaming + collision the sweep inverts directly), :math:`S` is the
within-group scattering gain (:ref:`pn-scattering` in
:doc:`discrete_ordinates`), :math:`B` is the reflective boundary
coupling (:ref:`bc-extraction`), and :math:`q` the fixed external /
fission source. Both :math:`S` and :math:`B` are **lagged gains** — the
sweep never re-scatters mid-sweep (cf. the variadic driver,
:ref:`bc-extraction-variadic-driver`).

The load-bearing observation is the **arity of the scattering gain**.
:math:`S\,\psi` depends on :math:`\psi` *only* through its
spherical-harmonic flux moments. Writing the moment-projection operator
:math:`M` (the :eq:`flux-moments` quadrature contraction, the SH frame's
analysis face ``frame.analysis`` from
:meth:`Quadrature.angular_frame(L)
<orpheus.numerics.quadrature.Quadrature.angular_frame>`)

.. math::
   :label: angular-windowing-moment-projection

   \phi_\ell^m(\vec r)
   \;=\; (M\psi)_\ell^m(\vec r)
   \;=\; \sum_{n=1}^{N} w_n \, Y_\ell^m(\hat\Omega_n)\,
         \psi_n(\vec r),
   \qquad 0 \le \ell \le L,\; |m| \le \ell,

the scattering source factors **through the moment boundary**:

* the **isotropic** :math:`\ell = 0` (P0) and the **(n,2n)** doubling
  terms (:ref:`pn-scattering`) need only the scalar flux
  :math:`\phi_0 \equiv \phi_0^0`;
* the **anisotropic** :math:`P_{\ell\ge 1}` term needs the higher
  moments :math:`\phi_\ell^m` up to the scattering order :math:`L`.

So the per-ordinate iterate :math:`\psi \in \mathbb{R}^N` is mapped, at
the very first thing the sweep's source assembly does, onto the
:math:`(L{+}1)(2L{+}1)`-dimensional moment space — and **nothing
downstream of that projection ever reads the discarded
:math:`N - (L{+}1)(2L{+}1)` angular degrees of freedom**. The iterate
carries strictly more than the iteration consumes. Angular windowing
holds the iterate at the consumed dimension: the persistent state is
the moment tensor

.. math::
   :label: angular-windowing-moment-iterate

   \phi \;\in\; \mathbb{R}^{(L+1)\times(2L+1)\times n_g \times n_x \times n_y}
   \quad(\texttt{HarmonicMomentFlux}),
   \qquad\text{not}\qquad
   \psi \;\in\; \mathbb{R}^{N \times n_g \times n_x \times n_y}
   \quad(\texttt{AngularFlux}).

.. vv-status: angular-windowing-moment-iterate documented

For :math:`N = 110` (Lebedev order 17) and :math:`L = 1`
(:math:`(L{+}1)(2L{+}1) = 6`) the angular dimension drops **18.3×**.


.. _sn-angular-windowing-factoring:

The scattering factoring — :math:`S_{\rm aniso} = \tfrac{1}{W}\,R\,\Lambda\,M`
------------------------------------------------------------------------------

Phase 5a's commit 1 (``93807aa``) makes the factoring *structural* so
that the windowed and full-angular paths share one source of truth. The
anisotropic in-scatter is the §9 operator composition (the §15.2
sum-of-tensor-products form,
:meth:`~orpheus.transport.operators.scattering.ScatteringOperator.build_aniso_source`)

.. (vv-status rationale) The R·Λ·M anisotropic-scattering factoring.
   A structural identity (associativity of the three-operator
   composition); the verifiable content is the bit-identity of the two
   evaluation arms, pinned by the de-risk probe Q2/Q3 (0 ULP) and the
   independent Bell & Glasstone hand reconstruction Q2b (1.5 ULP).
   Documented.
.. vv-status: angular-windowing-aniso-factoring documented

.. math::
   :label: angular-windowing-aniso-factoring

   Q^{\rm aniso}_n(\vec r)
   \;=\; \frac{1}{W}\,\bigl(R\,\Lambda\,M\,\psi\bigr)_n(\vec r),
   \qquad W = \sum_n w_n,

where (reading right to left):

* :math:`M` is the moment **projection** :eq:`angular-windowing-moment-projection`
  (the SH frame's analysis face ``frame.analysis``);
* :math:`\Lambda` is the per-:math:`\ell` block-diagonal scattering on
  moment space :math:`\Lambda = \sum_\ell P_\ell \otimes \Sigma_{s,\ell}`
  (:class:`~orpheus.transport.operators.scattering.LegendreMomentScattering`);
* :math:`R` is the addition-theorem **reconstruction** with the
  :math:`(2\ell+1)` factor
  (the SH frame's reconstruction face ``frame.reconstruction``);
* :math:`1/W` is the producer-side normalization applied at the
  :meth:`~orpheus.transport.operators.scattering.ScatteringOperator.apply` boundary.

The associativity :math:`(R\,\Lambda)\,M` is the whole trick. The
**full-angular** path evaluates the composition as written —
:math:`R\cdot\Lambda\cdot M(\psi)`: project, scatter, reconstruct,
consuming the §5.6 :attr:`kernel <orpheus.transport.operators.scattering.ScatteringOperator.kernel>`
``= frame.conjugate(Λ)`` as a single composed ``np.ndarray`` operator. The
**windowed** path's iterate bulk **is** the moments :math:`\phi = M\psi`,
so :math:`M` is *already done* and only the **moment → source** map
:math:`R\,\Lambda` remains: :math:`R\cdot\Lambda(\phi)`. The dispatch is on
the iterate type: the
:class:`~orpheus.transport.fields.angular_flux.AngularFlux` arm of
:meth:`ScatteringOperator.apply <orpheus.transport.operators.scattering.ScatteringOperator.apply>`
does the :math:`M` projection first; the
:class:`~orpheus.transport.fields.harmonic_moment_flux.HarmonicMomentFlux`
arm skips it.

.. note::

   **P4 — the windowed arm now takes the explicit typed edges.** Earlier
   the two arms shared the single ``np.ndarray`` primitive
   :meth:`~orpheus.transport.operators.scattering.ScatteringOperator._aniso_source_from_moment_values`
   (``frame.reconstruct_after(Λ)``). The Frame campaign's P4 phase
   re-expressed the windowed arm's :math:`R\,\Lambda` as the **explicit
   typed** carrier path
   ``frame.reconstruct(Λ.apply(φ)) / W`` — :math:`\Lambda` materialises a
   typed :class:`~orpheus.transport.source_sinks.harmonic_moment_source_sink.HarmonicMomentSourceSink`
   (the role-changing edge of the carrier grid,
   :ref:`scattering-carrier-grid`), which the frame's :math:`R` face then
   reconstructs to an :class:`AngularSourceSink`. It threads the *same*
   :math:`\Lambda` kernel and the *same* frame :math:`R` face as the fused
   path, so the two agree numerically; the ndarray
   ``reconstruct_after(Λ)`` primitive is **retained as the 0-ULP
   crosscheck oracle** the typed arm is pinned against, not deleted. See
   :ref:`scattering-carrier-grid` for the explicit-vs-fused design choice.

This factoring **retired** the per-:math:`\ell`
``_PerLegendreOrderScattering`` kernel, which recomputed :math:`M\psi`
independently for every Legendre order — an :math:`L`-fold redundant
projection (aggressive-retirement discipline).

.. admonition:: The :math:`Y_0^0 = 1` convention — the scalar flux is read off
   :class: note

   ORPHEUS uses **unnormalized real harmonics** (the
   "no-:math:`4\pi/(2\ell+1)`-prefactor" convention,
   :ref:`spherical-harmonics`), under which :math:`Y_0^0 = 1` *exactly*.
   Therefore the :math:`\ell = 0` moment **is** the scalar flux,

   .. math::

      \phi_0 \;=\; \sum_n w_n \, Y_0^0(\hat\Omega_n)\,\psi_n
                 \;=\; \sum_n w_n \, \psi_n
                 \;=\; (\texttt{integrate\_angular}\;\psi),

   read directly off the moment tensor's :math:`\ell{=}0` block with no
   rescale (``phi_moments.l_block(0)[0]``). This is what lets the
   windowed P0 + (n,2n) fast path consume the moments with **zero**
   conversion arithmetic, and what makes the eigenvalue outer's scalar
   flux bit-identical to the full-angular ``integrate_angular`` (the
   de-risk Q1 below proves :math:`\Phi[0,0] = \texttt{integrate\_angular}`
   bit-for-bit).


.. _sn-angular-windowing-geometry-restriction:

Why 2-D Cartesian only — the curvilinear seed obstruction
---------------------------------------------------------

Windowing is valid **only** where the within-group sweep is a *direct*
solve that does not seed from the previous iterate's per-ordinate
:math:`\psi`. There are three regimes, and only one admits it:

.. list-table:: Where angular windowing applies
   :header-rows: 1
   :widths: 26 16 58

   * - Path
     - Windowed?
     - Why
   * - **2-D Cartesian DD** (SI inner)
     - **yes**
     - The diamond-difference wavefront sweep is a direct forward
       substitution down the upwind DAG; it inverts :math:`L+C` from
       :math:`q + S\phi + B\psi_\partial` with **no** interior-iterate
       seed. The bulk seed (``initial_guess``) threads through harmlessly
       (the 2-D sweep ignores it). The moment iterate is sufficient.
   * - **Curvilinear 1-D** (sphere / cylinder)
     - **no**
     - The Morel–Montry **Carlson coupled-pole** closure seeds the
       :math:`\mu = -1` starting-direction angular flux from the
       *previous iterate's per-ordinate* :math:`\psi` (the curvilinear
       angular-redistribution recursion is initialized from the inward
       radial sweep's last iterate). A moment tensor cannot reconstruct
       that per-ordinate seed — windowing would lose the closure's
       starting data. Curvilinear stays full-angular.
   * - **Krylov** (any geometry)
     - **no**
     - GMRES iterates the full bulk vector :math:`\psi`; the Krylov
       subspace is built from full-angular matvecs. There is no moment
       sub-iterate to hold.

The gate is the genuine predicate ``sn_mesh.is_cartesian and
sn_mesh.ndim == 2`` (the C5.4 / #225 sharpening of the earlier
``reduced is None`` proxy, which was *also* true at 3-D Cartesian and
would have silently moment-windowed a 3-D solve — vv Mode 9): the
curvilinear meshes carry a non-``None`` ``reduced`` moment-reduction
descriptor and are excluded, and so is 3-D Cartesian (the in-sweep
moment emit is a 2-D kernel). The windowed product ``P @ A.inverse()``
(retyped in :ref:`windowing-retyped` below) — which the SI driver now
holds and applies **directly**, the ``_maybe_window`` factory returning
the plain sweep off that gate (#226 taxonomy step 3,
:ref:`inverse-application-driver`) — is therefore **never even
constructed** off the 2-D-Cartesian gate, so there is no illegal state
to mistype.

The restriction is **interior-bulk-only** in a second sense: the
interior-cochain biproduct :eq:`wavefront-cochain-biproduct`
:math:`C^1 = C^1_{\rm int} \oplus C^1_\partial` keeps the **boundary
trace** :math:`C^1_\partial` (the :class:`AngularBoundaryFlux` summand) a
distinct, *un-reduced* per-ordinate object. The reflective :math:`B`
coupling reads the full per-ordinate face trace via the typed
:math:`\iota_*` / :math:`\iota^*` exchange — windowing reduces only the
interior bulk and never touches the trace. The
``test_2d_windowed_si_reflective_trace_is_nonzero`` guard pins this (a
windowing that zeroed the trace would be a dropped reflective coupling,
invisible to the interior-only scalar-flux snapshot).


.. _sn-angular-windowing-bit-identity:

Bit-identity of the source, principled-equivalence of the convergence test
---------------------------------------------------------------------------

The carve has a clean correctness story split in two
(``vv-principles`` § "Bit-identity vs principled-equivalence").

**The per-step source is bit-identical.** A de-risk probe
(``derivations/diagnostics/diag_p5a_moment_consuming_scatter.py``) proved
the **moment arm** :math:`S(M\psi)` equals the **full-angular arm**
:math:`S(\psi)` **bit for bit** (``np.array_equal``, 0 ULP) before any
production code was written. Because the moment-projection operator
:math:`M` inside the windowed resolvent is built from the **same**
quadrature harmonics :math:`Y` and weights :math:`w` the scattering
operator uses internally, the stored moments equal :math:`S`'s own
internal projection of the same :math:`\psi` term-for-term — so the
per-sweep re-projection is not just *elided*, its result is *reproduced
exactly*. The probe was cross-checked against a **structurally
independent** Bell & Glasstone §1.6 hand reconstruction of the P1 source
(every factor — :math:`w_n`, :math:`Y`, the :math:`(2\ell+1)=3`, the
:math:`1/W` — written out by hand from the material data, *not* via the
project's projection primitives), which agreed at ~1.5 ULP (the expected
floating-point distance for an independent reduction order). This is the
L11 structural-independence guard the bit-exact comparison alone lacks.

.. list-table:: De-risk probe — moment arm vs full-angular arm
   :header-rows: 1
   :widths: 30 14 28 28

   * - Probe (config: 2-D P1 het, LS-S4)
     - Groups
     - Result
     - Verdict
   * - **Q1** :math:`\Phi[0,0] = \texttt{integrate\_angular}`
       (the :math:`Y_0^0 = 1` scalar-flux read)
     - 2G
     - max :math:`|\Delta| = 0`, ``np.array_equal`` True
     - bit-exact (0 ULP)
   * - **Q2** :math:`R\Lambda(\phi)` vs
       :math:`R\Lambda M(\psi)` (the aniso :math:`\ell\ge 1` arm)
     - 2G
     - max ULP :math:`= 0`
     - bit-exact (0 ULP)
   * - **Q3** full :math:`S(M\psi)` vs :math:`S(\psi)`
       (end-to-end source)
     - 2G **and** 4G
     - max ULP :math:`= 0` (both)
     - bit-exact (0 ULP)
   * - **Q2b** vs INDEPENDENT Bell & Glasstone hand
       reconstruction (L11 structural-independence ground)
     - 2G
     - max rel :math:`= 3.4\times10^{-16}`
     - principled (~1.5 ULP)
   * - **Q4** non-degeneracy: aniso :math:`> 0`,
       P1 :math:`\ne` P0, asymmetric :math:`\Sigma_s` (ERR-002 ground)
     - 2G
     - aniso max :math:`= 0.49`, :math:`|\Sigma_{s,0}-\Sigma_{s,0}^\top| = 0.1`–:math:`0.18`
     - non-degenerate (gate can see Mode 6)

**The convergence test moves to moment space (principled-equivalence).**
The *only* non-bit-identical change is the SI stopping criterion. The
full-angular path tested :math:`\lVert\psi_{k+1} - \psi_k\rVert_2`; the
windowed path necessarily tests :math:`\lVert\phi_{k+1} -
\phi_k\rVert_2`. This is the **physically-meaningful** SI criterion —
scattering iterates the *moments*, so the moment :math:`L^2` is the
natural convergence norm — i.e. the change is *more* principled, not a
regression. Because the stopping point shifts by a fraction of a
tolerance, the converged value agrees with the full-angular path within
``SAFETY × conv_tol``: the measured drift is **2-D eigenvalue
:math:`k_{\rm eff}` :math:`2.4\times10^{-14}` relative** and **scalar
flux :math:`6.3\times10^{-12}` relative**, both well inside the
regression framework's ``kind="iterative"`` gate (drift bounded by
:math:`(\text{iteration count})\times(\text{condition number})\times
\text{ULP}`, criterion 3 of the principled-equivalence test). The
eigenvalue **outer** converges on the scalar flux, which is bit-identical
via :math:`\phi_0` (the :math:`Y_0^0 = 1` read) — so *only* the inner SI
stopping point shifts; the outer power iteration is untouched.

The one *correctness* anchor (not merely *equivalence*) is the
closed-form homogeneous eigenvalue: the windowed 2-D eigenvalue solve
reproduces :math:`k_\infty = \nu\Sigma_f / \Sigma_a` (the closed-form
pillar — MMS does **not** prove eigenvalues), and the
``test_2d_p1_aniso_moment_path_carries_signal_and_si_krylov_agree`` gate
cross-checks the windowed SI flux against the **full-angular Krylov**
flux (a genuinely independent reference, since Krylov is never
windowed), confirming the windowing does not silently drop the
:math:`\ell\ge 1` moments — the central trap.


.. _sn-angular-windowing-honest-scope:

Honest scope — a persistent-iterate + typed-state win, NOT yet a peak win
-------------------------------------------------------------------------

.. warning::

   Phase 5a reduces the **persistent** iterate storage and makes the SI
   state honest. It does **NOT** by itself deliver the full peak-memory
   reduction. State this carefully — the 18.3× number describes the
   *held* iterate, not the *peak*.

   * **What 5a wins.** The held + warm-started ``_psi_typed`` carried
     across the entire eigenvalue solve (and the convergence-test copy)
     shrinks by :math:`N / (L{+}1)(2L{+}1)` — measured **18.3×** at
     :math:`N = 110`, :math:`L = 1`. The iterate **type** also becomes
     honest: the SI state *is* the moments
     (:class:`HarmonicMomentFlux`), so the representation no longer
     over-claims an angular resolution the iteration never uses.
   * **Why the peak win is modest.** The **per-sweep transient**
     full-angular machinery still dominates the peak: the resolvent's
     swept output, the ``psi_x`` / ``psi_y`` interior cochain
     :math:`C^1_{\rm int}` (the cochain section above), and the per-octant
     ``.copy()`` buffers — several full-angular-sized arrays that
     storage-A still materializes within a single sweep. Measured
     **peak** reduction is only **~1.2×** for a :math:`12\times8` config
     (the held :math:`2\times` iterate restored to full-angular size
     against the windowed ``tracemalloc`` peak;
     ``derivations/diagnostics/diag_p5a_peak_memory.py``).

   The per-sweep transient has two components, eliminated by two later
   phases. **Phase 5b** (interior-face storage-B — the rolling
   moving-frontier window over the wavefront anti-diagonals, which never
   materializes the whole interior cochain at once) cuts the *interior
   face* cochain transient and is the **3-D enabler** (a 3-D wavefront's
   full interior cochain is prohibitively large to materialize; the moving
   frontier is the only tractable representation). **Phase 5c**
   (:ref:`sn-angular-windowing-in-sweep-accumulation`) cuts the
   *full-angular output* transient by accumulating the harmonic moments
   per anti-diagonal inside the sweep — the **linear peak win** (measured
   3.06×), and the realization of the full peak reduction this honest-scope
   warning deferred.

So Phase 5a is precisely: the **persistent-iterate** shrink + the
**typed-state** correctness win + the **foundation** for 5b/5c. Phase 5b
eliminates the **interior-face transient** and is the **3-D** enabler;
Phase 5c eliminates the **full-angular output transient** (the linear peak
win). The three together deliver the asymptotic peak reduction the moment
iterate makes possible; 5a alone delivers the type and the
persistent-storage win.


.. _sn-angular-windowing-implementation:

Implementation map
-------------------

* The solver's :func:`_maybe_window <orpheus.sn.solver._maybe_window>`
  builds the typed windowed composition ``P @ A.inverse()`` — the
  :class:`~orpheus.sn.operators.windowing.WindowedSweep` — over the
  within-group forward :math:`A`'s inverse (the Jacobi
  :math:`(L+C)^{-1}`, or the reified splitting matrix :math:`M^{-1}`,
  :class:`~orpheus.sn.operators.scheduled_invertible.ScheduledInvertibleOperator`
  — see :ref:`si-gauss-seidel-reification` in :doc:`discrete_ordinates`)
  and hands it to the :class:`~orpheus.numerics.iteration.SourceIteration`
  driver, which **applies** it directly (#226 taxonomy step 3 —
  :ref:`inverse-application-driver`; the transitional
  ``_MomentWindowedResolvent`` ``.solve`` adapter that once wrapped it is
  gone). Its analysis factor :math:`P` is sourced from the scattering
  operator's **own** frame, which is what makes the stored moments equal
  :math:`S`'s internal projection term-for-term. The composite reports
  :attr:`~orpheus.numerics.operator.LinearOperator.is_invertible`
  ``= False`` — no round-trip promise (its :math:`P`
  factor is a coisometry) — and accepts the driver's ``initial_guess``
  kwarg accepted-and-ignored (:meth:`WindowedSweep.apply
  <orpheus.sn.operators.windowing.WindowedSweep.apply>`; the multi-D walk
  has no bulk-seed consumer). Gated on ``sn_mesh.is_cartesian and
  sn_mesh.ndim == 2``.
* The **eigenvalue** inner
  (:meth:`SNSolver._solve_source_iteration <orpheus.sn.solver.SNSolver._solve_source_iteration>`)
  is reconstruction-free: it returns the scalar flux read off the
  :math:`\ell{=}0` moment block (``psi_typed.bulk.l_block(0)[0]``, the
  :math:`Y_0^0 = 1` identity), since the outer power iteration rebuilds
  from the scalar flux.
* The **fixed-source** inner
  (:func:`~orpheus.sn.solver.solve_sn_fixed_source`'s SI path) adds a
  **one-shot** full-angular reconstruction for ``Solution.angular_flux``
  — it re-evaluates the converged source ``q + Σ gains·ψ`` once through
  the *un-wrapped* base resolvent (``base_resolvent.solve(...)``).
  Fixed-source must return the full per-ordinate field, and because
  :math:`S`/:math:`B` consume the moments *equal* the full angular's
  moments (de-risk proven), the reconstructed source is the same and one
  sweep reproduces the converged iterate by the fixed point — so the
  reconstruction is bit-identical to the un-windowed converged
  :math:`\psi`.

Cross-references: the moment math :eq:`flux-moments` / :eq:`pn-scatter`
and the :math:`Y_\ell^m` convention live in :ref:`pn-scattering`
(:doc:`discrete_ordinates`); the projection :math:`M`
(equation :math:`\phi_\ell^m = \sum_n w_n Y_\ell^m \psi_n`) and the
addition-theorem reconstruction :math:`R` are the spherical-harmonic
:class:`~orpheus.numerics.frame.GalerkinFrame`'s ``analysis`` /
``reconstruction`` faces (:ref:`galerkin-projection`); the interior-face
cochain Phase 5a is orthogonal to is :ref:`wavefront-flux-cochain`; the
SI fixed point it reorganizes is the within-group inner of
:ref:`eigenvalue-posing`.


.. _sn-angular-windowing-in-sweep-accumulation:

Per-anti-diagonal moment accumulation — dropping the full-angular transient (Phase 5c)
--------------------------------------------------------------------------------------

Wave O step #205 **Phase 5c** (commit ``c7be111``, 2026-06-07) closes the
gap Phase 5a's :ref:`honest-scope warning
<sn-angular-windowing-honest-scope>` left open: the **per-sweep
full-angular transient**. Phase 5a held the *persistent* iterate as
moments (the 18.3× shrink) but still **materialized the full
per-ordinate angular field** :math:`\psi \in (N, n_g, n_x, n_y)` inside
*every* sweep, then projected it to moments **post-hoc** — a flat reduce
:eq:`angular-windowing-moment-projection` applied once at the end of the
sweep by the SH frame's analysis face ``frame.analysis.apply`` (the then
Phase-5a moment-windowed resolvent's ``solve`` called ``base.solve`` then
``self._projection.apply(full.bulk.values)``).
That transient full-angular array is the dominant peak-memory cost the
5a warning named.

Phase 5c moves the projection **into** the windowed anti-diagonal walk:
each topological level accumulates the harmonic moment tensor *directly*,
octant-by-octant into one shared global buffer, so the full angular
field is **never materialized** in the windowed iterate. The
``base.solve`` → flat-``apply`` post-projection **leaves production**; the
fused moment emit is what survives — reached at Phase 5c through the
resolvent's ``solve``, at #226 step 2 through the base resolvent's
``solve_moments``, and since #226 step 3 through the typed product's
:meth:`WindowedSweep.apply <orpheus.sn.operators.windowing.WindowedSweep.apply>`
that the :class:`~orpheus.numerics.iteration.SourceIteration` driver
applies directly (:ref:`windowing-retyped` and
:ref:`inverse-application-driver`, below).

.. admonition:: Key Facts (in-sweep moment accumulation)
   :class: tip

   - **The projection moves into the sweep.** Where 5a swept the full
     :math:`\psi` then flat-reduced :math:`\phi_\ell^m = \sum_n w_n
     Y_\ell^m \psi_n` once post-hoc, 5c accumulates the moments per
     anti-diagonal *during* the walk:
     :math:`\phi_\ell^m[\text{cells}_k] \mathrel{+}= \sum_{n\in
     \text{octant}} w_n Y_\ell^m(\hat\Omega_n)\,\psi_n[\text{cells}_k]`.
     The full :math:`(N, n_g, n_x, n_y)` field is never allocated for
     the windowed iterate.
   - **Measured peak-memory win: 3.06×.** On S8 / 4g / :math:`24\times24`
     the windowed ``solve`` peak drops 2.26 MB → 0.74 MB — the 1.47 MB
     full-angular transient is eliminated (moment tensor 0.111 MB). The
     win **grows with angular order**: the eliminated transient is
     :math:`N\,n_g\,n_x\,n_y`; the moment tensor is fixed at
     :math:`(L{+}1)(2L{+}1)\,n_g\,n_x\,n_y`.
   - **Principled-equivalence, NOT bit-identity.** The cross-octant
     ``+=`` reorders the ordinate sum vs the flat single-reduce; IEEE-754
     addition is non-associative, so the moments drift at ULP level.
     The per-cell :math:`w\,Y\,\psi` fold is **term-for-term identical**
     to the SH frame's analysis face ``frame.analysis.apply`` — only the
     accumulation *order* differs. Measured max-relative drift
     :math:`2.74\times10^{-16}` (:math:`\le 4N\varepsilon`, 4 ULP).
   - **The scalar is subsumed.** The scalar flux **is**
     :math:`\phi_0^0 = \texttt{moments}[0,0]` (:math:`Y_0^0 = 1`), read
     off the moment tensor — there is no separate scalar reduction
     (``coding-elegance`` Pattern 2: one source of truth).
   - **The fuller view is retained as a verification oracle.** The pre-5c
     "full-angular solve + flat ``frame.analysis.apply``" path is kept
     reachable and pins the optimized path
     (``feedback_aggressive_retirement`` — the "verification oracle"
     exception to retirement).
   - **2-D Cartesian only; the rest is untouched.** Gated on the genuine
     ``sn_mesh.is_cartesian and sn_mesh.ndim == 2`` (the C5.4 / #225
     sharpening of 5a's ``reduced is None`` proxy — the proxy was also
     true at 3-D Cartesian); ``transport_sweep`` *raises*
     on a 1-D mesh given a moment projection (illegal-states-
     unrepresentable). Curvilinear / Krylov stay full-angular; both
     one-shot ``Solution.angular_flux`` reconstructions stay separate
     full sweeps.


.. _sn-angular-windowing-in-sweep-transformation:

The transformation — post-hoc reduce → in-sweep accumulate
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Phase 5a's resolvent produced its moment iterate in **two arithmetic
stages**: (1) the wavefront walk materialized the full per-ordinate field
:math:`\psi \in (N, n_g, n_x, n_y)`, writing each anti-diagonal's
cell-average into the global angular buffer
(``angular_flux_octant[:, :, ii, jj] = psi_avg``); then (2) a flat
post-sweep reduce collapsed *all* :math:`N` ordinates at once,

.. math::

   \phi_\ell^m(\vec r)
   \;=\;\sum_{n=1}^{N} w_n\, Y_\ell^m(\hat\Omega_n)\,\psi_n(\vec r)
   \qquad\bigl(\texttt{frame.analysis.apply}, \;
   \texttt{einsum}\;\texttt{"n,nlm,n...->lm..."}\bigr),

which is :eq:`angular-windowing-moment-projection` again, evaluated once.
Stage (1)'s full-angular array is the transient that dominates the peak.

Phase 5c fuses the two stages. The :class:`_MovingFrontier
<orpheus.sn.loss_representation.sweep_graph._MovingFrontier>` walk
(:meth:`SweepDependencyGraph.walk_windowed
<orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.walk_windowed>` — at 5c the
solve-direction ``apply_windowed``, collapsed into the level-op walk at
S6.4(e)) already
visited every cell exactly once, anti-diagonal by anti-diagonal, with the
per-level cell-average :math:`\psi_n[\text{cells}_k]` in hand. Phase 5a
threw that local quantity into a global angular buffer and re-read it
later; 5c **projects it on the spot** and discards it. For each level
:math:`k` (anti-diagonal ``cells_k = (ii, jj)``) of each in-plane octant,

.. (vv-status rationale) The in-sweep per-anti-diagonal moment
   accumulation. The verifiable claim is the fuller-view-oracle
   equivalence: the in-sweep accumulation equals the flat post-sweep
   frame.analysis.apply of the same swept psi within the
   reduction-order drift bound (criterion 3 of bit-identity-vs-
   principled-equivalence). Pinned by
   test_2d_windowed_product_equals_post_projection (renamed from
   test_2d_windowed_moments_in_sweep_equal_post_projection when the
   solve_moments cross-reach retired into the typed product, #226 step 2),
   anchored to the structurally-independent SI≡Krylov-full scalar
   cross-check.
.. vv-status: harmonic-moment-projection implemented

.. math::
   :label: harmonic-moment-projection

   \phi_\ell^m[\text{cells}_k]
   \;\mathrel{+}=\;
   \sum_{n\in\text{octant}} w_n\, Y_\ell^m(\hat\Omega_n)\,
       \psi_n[\text{cells}_k],
   \qquad
   \texttt{moment\_buf[:, :, :, ii, jj]}
   \mathrel{+}=
   \texttt{einsum("nlm,ngd,n->lmgd",}\,
       Y_{\rm oct},\,\psi_{\rm avg},\,w_{\rm oct}\texttt{)},

accumulated into **one shared global** ``moment_buf`` of shape
:math:`(L{+}1, 2L{+}1, n_g, n_x, n_y)`. The four in-plane octants are
swept sequentially (the octant loop — since S6.4(b) the shared
``_OctantWalk.sweep_group`` frame, then ``sweep_octant_group`` — each octant
completing its full frontier walk), so the global buffer receives a **cross-octant
``+=``**: octant 1's complete contribution, then octant 2's, and so on.
The output branch (since S6.4(e) inside the ``_CellSolve`` level
operation; at 5c inside ``apply_windowed``) is a single
two-way switch at the per-level output site — *angular mode* writes
``angular_flux_octant[:, :, ii, jj] = psi_avg`` and accumulates the scalar
(``scalar_flux_buf``); *moment mode* does the :eq:`harmonic-moment-projection`
``+=`` instead. The frontier walk and the cell kernel are **untouched**,
which is exactly why the ``window ≡ full-field`` angular-mode oracle
(:ref:`Phase 5b <wavefront-flux-cochain>`) stays bit-identical: 5c adds a
branch only at the *consumer* of ``psi_avg``, never in its *production*.

The reduction-order contrast is the whole story:

.. list-table:: Post-hoc reduce (5a) vs in-sweep accumulate (5c)
   :header-rows: 1
   :widths: 22 39 39

   * -
     - **Phase 5a — post-hoc flat reduce**
     - **Phase 5c — in-sweep accumulation**
   * - Full :math:`\psi` materialized?
     - **Yes** — written to the global angular buffer, then re-read
     - **No** — projected per level and discarded
   * - Reduction shape
     - one flat ``einsum`` over all :math:`N` at once
       (``"n,nlm,n...->lm..."``)
     - per-level ``+=`` per octant
       (``"nlm,ngd,n->lmgd"``), cross-octant accumulate
   * - Per-cell arithmetic
     - :math:`\sum_n w_n Y_\ell^m \psi_n`
     - :math:`\sum_n w_n Y_\ell^m \psi_n` (**term-for-term identical**)
   * - FP reduction tree
     - one fixed order
     - reordered by the octant grouping ⟹ ULP drift
   * - Peak working set
     - moment iterate **+** full-angular transient
     - moment iterate **only**


.. _sn-angular-windowing-in-sweep-equivalence:

Why it is principled-equivalence, not bit-identity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Phase 5a's per-step source was **bit-identical** to the full-angular arm
(0 ULP — :ref:`sn-angular-windowing-bit-identity`): the moments *stored*
equalled :math:`S`'s *own* internal projection of the same :math:`\psi`,
because the flat reduce and :math:`S`'s internal reduce shared the same
single reduction order. Phase 5c **necessarily** breaks that bit-identity
— and does so for a clean, documented reason rooted in the
``vv-principles`` § "Bit-identity vs principled-equivalence" three-criteria
test. The change is **accepted** because all three criteria hold:

#. **Principled at every step — the intermediate is named.** Each
   per-anti-diagonal partial moment :math:`\sum_{n\in\text{octant}} w_n
   Y_\ell^m(\hat\Omega_n)\,\psi_n[\text{cells}_k]` is **the octant's
   contribution to the harmonic moment** :math:`\phi_\ell^m` — a
   reactor-physics quantity (the per-direction-class flux moment), not a
   "whatever the reduction order produced" array. The cross-octant ``+=``
   composes named partial moments. This is exactly the
   unnamed-intermediate → named-intermediate move the criterion blesses
   (cf. the issue-#169 ``compute_keff`` worked example: per-group
   production rate is principled; the flat cell-product array is not).
#. **Verified against a structurally-independent reference — not
   old-vs-new ULP alone.** Old-vs-new distance is *necessary but not
   sufficient* (both arms could share a systematic offset). The
   structurally-independent ground is the **full-angular Krylov** scalar:
   ``test_2d_p1_aniso_moment_path_carries_signal_and_si_krylov_agree``
   cross-checks the windowed-SI moment :math:`\ell{=}0`
   (``scalar_flux()`` = :math:`\phi_0^0`) against the full-angular Krylov
   flux (Krylov is **never** windowed — it iterates the full bulk vector),
   and the closed-form homogeneous eigenvalue :math:`k_\infty =
   \nu\Sigma_f/\Sigma_a` (the closed-form pillar; MMS does **not** prove
   eigenvalues) anchors the eigenvalue reuse.
#. **The drift is FP-non-associativity, dimensionally explainable.** For
   a single-step computation the bound is :math:`(\text{reduction
   depth})\times\varepsilon`. The reduction depth is the ordinate count
   :math:`N`; the cross-octant ``+=`` reorders the **same** :math:`N`
   summands. The de-risk probe (``derivations/diagnostics/diag_p5c_moment_accum.py``,
   git-excluded, promoted into the permanent gate then deleted per the
   diagnostic-promotion policy) measured max-relative drift
   :math:`2.74\times10^{-16}` — about :math:`78\times` under the
   :math:`4N\varepsilon \approx 2.1\times10^{-14}` bound (:math:`N = 24`
   ordinates, 4× headroom for the partial-sum nesting within octants). A
   drift *above* :math:`4N\varepsilon` would signal an algorithmic change
   masquerading as FP noise — a wrong octant-:math:`Y` slice (Mode 2), a
   missing weight (Mode 3), or an :math:`\ell/m` index drift (Mode 5) —
   and the gate would catch it.

The scalar flux is **subsumed**, not separately reduced. Because
:math:`Y_0^0 = 1` (ORPHEUS unnormalized-harmonic convention,
:ref:`spherical-harmonics`), the :math:`\ell{=}0` slot **is** the scalar
flux, :math:`\phi_0^0 = \sum_n w_n \psi_n` (the :math:`Y_0^0 = 1` read of
:ref:`sn-angular-windowing-factoring`). The moment-mode sweep therefore
returns ``(moment_buf, None)`` — no separate ``scalar_flux_buf``
accumulation — and the eigenvalue inner reads
``moment_buf[0, 0]`` for the outer power iteration's scalar flux. The
windowed scalar is identical to the angular-mode scalar up to the same
reduction-order drift; the existing ``scalar_flux_buf`` einsum
``"ngd,n->gd"`` is literally the :math:`\ell{=}0` case of the moment
einsum ``"nlm,ngd,n->lmgd"`` with :math:`Y_0^0 = 1`.


.. _sn-angular-windowing-in-sweep-oracle:

The fuller-view verification oracle
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Phase 5c **relinquishes a fuller view** — the materialized full-angular
field that the post-hoc reduce consumed. By the
``feedback_aggressive_retirement`` "verification oracle" exception, an
optimization that gives up a fuller view of a concept **keeps that fuller
view reachable** as a verification oracle that pins the optimized path.
The pre-5c "full-angular ``base.solve`` then flat
``frame.analysis.apply``" path is exactly
that fuller view: it is *not* deleted, and — since #226 step 2 — it is
the **inherited**
:meth:`OperatorProduct.apply <orpheus.numerics.operator.OperatorProduct.apply>`
body of the windowed composition itself (``P.apply(A⁻¹.apply(rhs))``, the
un-fused factor-by-factor evaluation the
:class:`~orpheus.sn.operators.windowing.WindowedSweep` overrides). The
permanent gate ``test_2d_windowed_product_equals_post_projection``
asserts the fused in-sweep
:meth:`WindowedSweep.apply <orpheus.sn.operators.windowing.WindowedSweep.apply>`
result equals that deforested oracle of the same swept
:math:`\psi`, within the de-risk bound, over the **full moment tensor —
including the** :math:`\ell\ge 1` **block**.

The :math:`\ell\ge 1` coverage is the load-bearing reason this oracle
exists. The :math:`\ell{=}0` scalar cross-check (the SI ≡ Krylov-full
test, criterion 2 above) is **blind to** :math:`\ell\ge 1`: a windowing
that silently swapped or dropped a higher moment would converge the right
scalar flux while corrupting the anisotropic source. The oracle pins the
:math:`\ell\ge 1` block that the scalar cross-check cannot see — it is the
**Mode-5** (:math:`\ell/m` index drift) and **Mode-2** (wrong octant-:math:`Y`
slice) catcher.

.. warning::

   The oracle and the system-under-test share the **same** cell kernel
   (:meth:`~orpheus.sn.scheme.DiamondDifference.cell_kernel_batch`)
   and the **same** :math:`Y` / ``weights``, differing only in reduction
   *order*. They are therefore **procedurally** independent, **not
   structurally** independent (``vv-principles`` § "Structural
   independence"). The oracle gate alone is a same-math-rearranged
   comparison and is **not sufficient** on its own — it MUST be paired
   with the structurally-independent SI ≡ Krylov-full scalar anchor
   (criterion 2) and the closed-form :math:`k_\infty` eigenvalue. The
   oracle's *value-add* is precisely the :math:`\ell\ge 1` block the
   structurally-independent scalar anchor is blind to; its
   anti-degeneracy leg asserts :math:`\max|\phi_{\ell\ge1}| > 10^{-3}\,
   \max|\phi_0|` so the pinned higher-moment block is genuinely non-zero
   on the canonical config.

The metric the gate uses is the **scale-relative** drift
:math:`\max|\phi^{\rm SUT} - \phi^{\rm oracle}| / \max|\phi^{\rm oracle}|
\le 4N\varepsilon`, **not** an element-wise
``assert_array_almost_equal_nulp``. The moment tensor spans ~3 orders of
magnitude (the :math:`\ell{=}0` scalar dominates; the :math:`\ell\ge 1`
blocks are :math:`\sim 10^{-3}\times`), so an element-wise ULP comparison
would inflate a machine-:math:`\varepsilon` *absolute* difference on a
small :math:`\ell\ge 1` element into hundreds of ULP even when the reorder
is pure FP noise. The scale-relative drift is the principled-equivalence
quantity (criterion 3).


.. _sn-angular-windowing-in-sweep-numerical-evidence:

Numerical evidence
~~~~~~~~~~~~~~~~~~~

.. list-table:: Phase 5c gates (commit ``c7be111``)
   :header-rows: 1
   :widths: 34 18 48

   * - Gate
     - Result
     - What it pins
   * - **De-risk drift** (canonical config
       ``2d_2g_p1_aniso_dd_8x4_het_si``: 8×4 het, vacuum-x /
       reflective-y, mixture B :math:`\bar\mu = 0.6` P1, S4 :math:`N = 24`)
     - **2.74e-16** rel
     - in-sweep accumulation ≡ flat post-projection, :math:`\le 4N\varepsilon`
       (4 ULP) — criterion 3 (FP-reorder, not a bug)
   * - **Peak memory** (S8 / 4g / :math:`24\times24`,
       ``diag_p5c_peak_memory`` tracemalloc)
     - **3.06×** (2.26 MB → 0.74 MB)
     - the 1.47 MB full-angular transient eliminated; moment tensor 0.111 MB
   * - **Windowing L1**
       (``test_2d_anisotropic_windowing.py``, incl. the new equivalence gate)
     - **4 ✓**
     - P1 ≠ P0 + SI ≡ Krylov-full; full-recon self-consistency
       (:math:`\Sigma w\psi = \phi`); reflective trace ≠ 0; in-sweep ≡ oracle
   * - **Eigenvalue 2-D** (``test_keff_2d.py``)
     - **19 ✓**
     - closed-form :math:`k_\infty`, P1-changes-:math:`k_{\rm eff}`,
       SI ≡ Krylov, Jacobi ≡ G-S (Mode 9), refinement convergence
   * - **Regression snapshot**
       (``test_dd_regression.py``, ``2d_2g_p1_aniso_dd_8x4_het_si``)
     - within ``SAFETY × conv_tol``
     - scalar flux drift 6920 ULP / :math:`9.81\times10^{-13}` rel
       (~10× headroom under the :math:`1.0\times10^{-11}` gate)
   * - **5b oracles + matvec consumer**
       (``test_2d_full_field_oracle.py``, window≡full graph,
       ``TestAnisoMomentSourcePath``, A2D-1 hash pin)
     - bit-identical / intact
     - the angular-mode oracle and the moment *consumer* (``S.apply``)
       are untouched — 5c changes only moment *production*

The **peak-memory win grows with angular order** because the eliminated
transient scales as :math:`N\,n_g\,n_x\,n_y` while the moment tensor is
fixed at :math:`(L{+}1)(2L{+}1)\,n_g\,n_x\,n_y`: at S4 :math:`N = 24`,
:math:`L = 1` the moment tensor is :math:`6/24 = 1/4` the angular size;
at S8 :math:`N = 80` it is :math:`6/80 \approx 1/13`; a Lebedev-order-17
:math:`N = 110` it is :math:`6/110 \approx 1/18`. The transient is the
**linear** peak win Phase 5a's honest-scope warning deferred to "Phase
5b/5c": 5a delivered the persistent-iterate shrink, 5b (the
:ref:`moving-frontier interior cochain <wavefront-flux-cochain>`) cut the
interior face transient, and 5c cuts the full-angular *output* transient.


.. _sn-angular-windowing-in-sweep-honest-scope:

Honest scope — what 5c does and does NOT do
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. warning::

   Phase 5c eliminates the **per-sweep full-angular transient** of the
   *windowed SI* path. State the boundaries carefully:

   * **The persistent iterate was already moments (5a).** 5c does not
     shrink the held iterate — that was 5a's 18.3×. 5c shrinks the
     *transient* the held-moment iterate's sweep still materialized.
   * **Both** ``Solution.angular_flux`` **reconstructions stay full
     sweeps.** The eigenvalue final reconstruction and the fixed-source
     windowed one-shot reconstruction each re-run a *separate*
     full-angular ``transport_sweep`` to return the user-facing
     :math:`(N, n_g, n_x, n_y)` field. They are **untouched** — the
     user-facing angular flux is bit-identical (the
     :math:`\Sigma w\psi = \phi` self-consistency gate and the step-3
     ``np.array_equal`` reconstruction pin both hold).
   * **Krylov, 1-D, and curvilinear stay full-angular.** Krylov iterates
     the full bulk vector (no moment sub-iterate); the curvilinear
     Morel–Montry Carlson coupled-pole seed reads the per-ordinate
     iterate at :math:`\mu = -1` (lesson L21), which the moment tensor
     cannot carry. Both are gated out by the genuine
     ``sn_mesh.is_cartesian and sn_mesh.ndim == 2`` (C5.4 / #225 — 3-D
     Cartesian is excluded too); ``transport_sweep`` *raises* if a moment
     projection reaches a 1-D mesh, so the unwindowable regime is
     unrepresentable, not merely unreached.

So the three-phase arc is complete: **5a** holds the persistent iterate
as moments (typed-state + persistent-storage win); **5b** carries the
interior face cochain on a rolling moving frontier (interior-transient
elimination + 3-D enabler); **5c** projects the swept flux to moments
in-sweep (full-angular *output*-transient elimination — the linear peak
win). Together they realize the asymptotic peak reduction the moment
iterate makes possible.


.. _sn-angular-windowing-in-sweep-implementation:

Implementation map (Phase 5c)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The moment OUTPUT mode threads as an **optional projection**
(``moment_frame=None``), mirroring the ``reflect=None``
dependency-injection idiom of :func:`_sweep_scheduled <orpheus.sn.loss_representation._sweep_scheduled>` — **not** a boolean flag.

.. note::

   **Retyping (#226 step 2).** At Phase 5c the two output modes were
   **named methods** on the resolvent surface — ``solve`` (full angular)
   vs ``solve_moments`` (moments). That public ``solve_moments`` — a
   method whose output-mode argument silently changed the operator's
   *codomain* — was a composition wearing a config, and it was retired.
   The moment emit is now the typed composition ``P @ A.inverse()``
   (:ref:`windowing-retyped`); the ``moment_frame`` kwarg survives only
   on the **private** ``_solve_timed_full_field`` body, the single
   application-context entry that the fused
   :meth:`WindowedSweep.apply <orpheus.sn.operators.windowing.WindowedSweep.apply>`
   drives. The map below is updated to that state.

* the solve-direction windowed walk (at 5c ``apply_windowed``; since
  S6.4(e) :meth:`SweepDependencyGraph.walk_windowed
  <orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.walk_windowed>` × the
  ``_CellSolve`` level operation) — the single two-way branch at the
  per-level output site (angular write vs the
  :eq:`harmonic-moment-projection` ``moment_buf`` ``+=``). The
  frontier walk and cell kernel are untouched.  The apply-direction walk
  (at 5c ``residual_windowed``; now ``walk_windowed`` × ``_CellResidual``
  — the Krylov matvec) is untouched — Krylov stays full-angular.
* the per-group octant frame (``sweep_octant_group`` then, since S6.4(b),
  ``_OctantWalk.sweep_group``) /
  :func:`_sweep_scheduled <orpheus.sn.loss_representation._sweep_scheduled>` /
  :func:`_sweep_jacobi <orpheus.sn.loss_representation._sweep_jacobi>` /
  :func:`transport_sweep <orpheus.sn.loss_representation.transport_sweep>` — thread the
  optional ``moment_projection`` (2-D Cartesian only; ``transport_sweep``
  raises on a 1-D mesh). Moment mode skips the per-octant angular
  allocation and ``_sweep_2d_scheduled`` returns ``(moment_buf, None)``.
* the private
  :meth:`InvertibleOperator._solve_timed_full_field
  <orpheus.sn.operators.streaming.InvertibleOperator._solve_timed_full_field>`
  body (duck-shared by
  :meth:`ScheduledInvertibleOperator._solve_timed_full_field <orpheus.sn.operators.scheduled_invertible.ScheduledInvertibleOperator._solve_timed_full_field>`)
  is the **single** application-context entry: given a ``moment_frame`` it
  emits harmonic moments, else the full angular field. ONE
  representation-sweep call per body for both output modes (since S6.5 the
  operator's own ``loss_representation.sweep`` — the same instance the
  matvec consumes); only the bulk **wrap** differs
  (:class:`~orpheus.transport.fields.angular_flux.AngularFlux` vs
  :class:`~orpheus.transport.fields.harmonic_moment_flux.HarmonicMomentFlux`).
  The former public ``solve_moments`` cross-reach — on
  ``InvertibleOperator`` **and** the dissolved ``_GaussSeidelResolvent``,
  plus the solver-side ``_solve_scheduled`` plumbing twin — retired into
  this one private body (#226 step 2, :ref:`windowing-retyped`).
* the fused entry is
  :meth:`WindowedSweep.apply <orpheus.sn.operators.windowing.WindowedSweep.apply>`
  (``P @ A.inverse()``), which calls that private body with its
  ``moment_frame``; the analysis factor
  :class:`~orpheus.sn.operators.windowing.BulkAnalysisOperator` carries
  the SH frame (whose basis
  :class:`~orpheus.numerics.basis.SphericalHarmonicBasis` ``L`` carries
  the truncation order). Since #226 step 3 the
  :class:`~orpheus.numerics.iteration.SourceIteration` driver **applies**
  that fused ``apply`` directly (:ref:`inverse-application-driver`); the
  transitional ``_MomentWindowedResolvent`` ``.solve`` adapter that once
  forwarded to it is gone.

Cross-references: the post-hoc reduce 5c replaces is
:eq:`angular-windowing-moment-projection`; the :math:`Y_0^0 = 1`
scalar-flux read is :ref:`sn-angular-windowing-factoring`; the
moving-frontier interior cochain 5c rides is :ref:`wavefront-flux-cochain`
(Phase 5b); the persistent-iterate shrink 5c completes is
:ref:`sn-angular-windowing-honest-scope` (Phase 5a). The
principled-equivalence framework is ``vv-principles`` § "Bit-identity vs
principled-equivalence".


.. _windowing-retyped:

Windowing retyped — the moment emit as a composition (#226 step 2)
------------------------------------------------------------------

Phases 5a–5c built the moment-windowed path as a pair of **named
methods** on the resolvent surface: ``solve`` returned the full angular
field, ``solve_moments`` returned the harmonic moments. The taxonomy
step-2 carve (#226) retired ``solve_moments`` — because a *public method
whose output-mode argument silently changes the operator's codomain is a
composition wearing a config*, not a mode of one morphism.

The two output modes do not share a codomain: ``solve`` lands in the
full-angular composite carrier :math:`V = V_{\rm bulk} \oplus V_\partial`,
while ``solve_moments`` landed in the *moment-bulk* composite
:math:`V^{\rm mom} = \Phi_{\rm bulk} \oplus V_\partial`. A change of
codomain is, by the two-layer law of :ref:`operator-algebra` ("two
operators, one substrate — never two views, one operator"), a
**different arrow**; and an arrow whose target differs from what
:math:`A^{-1}` produces is a *composition* of :math:`A^{-1}` with a
second arrow, never a configuration flag on :math:`A^{-1}`. The honest
object is therefore

.. math::

   \text{windowed} \;=\; P \circ A^{-1},
   \qquad
   P \;=\; \underbrace{M_{\rm frame}}_{\text{analysis on the bulk}}
           \;\oplus\;
           \underbrace{\mathrm{Id}}_{\text{on the trace}},

where :math:`A^{-1}` is the within-group forward's inverse
(:class:`~orpheus.sn.operators.sweep_operator.SweepOperator` on
:math:`L+C`, or on the reified splitting matrix :math:`M` —
:ref:`si-gauss-seidel-reification`) and :math:`M_{\rm frame}` is the
scattering frame's :attr:`~orpheus.numerics.frame.FrameBase.analysis`
face, the angular→moment reduction :math:`\phi_\ell^m = \sum_n w_n
Y_\ell^m \psi_n` (:eq:`harmonic-moment-projection`). The boundary trace
passes through **un-reduced**: windowing is interior-bulk-only (the
reflective :math:`B` coupling reads the full per-ordinate face trace —
:ref:`sn-angular-windowing-geometry-restriction`), so :math:`P` is the
identity on the :math:`V_\partial` summand.

The coisometry factoring of the windowed contract
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:math:`P` — the :class:`~orpheus.sn.operators.windowing.BulkAnalysisOperator`
— is a **block coisometry**, not an isomorphism. Its bulk face satisfies

.. math::

   \text{analysis} \circ \text{reconstruction} \;=\; 4\pi\,\mathrm{I}

under the no-prefactor SH convention (:ref:`spherical-harmonics`) — the
addition-theorem tight-frame identity, pinned by
``test_pi_R_is_4pi_identity_through_the_frame``. It is emphatically **not**
:math:`\mathrm{I}`: asserting :math:`\Pi R = \mathrm{I}` was the ERR-051
mistake — the coisometry carries the :math:`4\pi` frame constant, and a
test that hard-coded :math:`= \mathrm{I}` verified the *wrong* invariant.
In the other order :math:`\text{reconstruction} \circ \text{analysis}
\neq \mathrm{I}` (moments discard the ordinate-resolved angular content
above the truncation order :math:`L`), so :math:`P` is **not
invertible**. By the product's invertibility-closure law the composite
:math:`P \circ A^{-1}` therefore honestly reports ``is_invertible =
False`` (its ``P`` factor is structurally non-invertible), and makes *no
round-trip promise* — which is the whole point: the old ``solve_moments``
name suggested a *solve* (an inverse) where the object is a **projection
composed with an inverse**.

Fusion is an evaluation strategy, not a new operator
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The composition itself is
:class:`~orpheus.sn.operators.windowing.WindowedSweep` (an
:class:`~orpheus.numerics.operator.OperatorProduct` subclass), spelled by
the operator algebra ``P @ A.inverse()`` — the
:meth:`BulkAnalysisOperator.__matmul__ <orpheus.sn.operators.windowing.BulkAnalysisOperator.__matmul__>`
dispatch recognises a :class:`SweepOperator` right factor and returns the
fused product (mirroring the ``L + C`` fusion precedent, #261: the
dispatch is one-directional on the specific leaf).
:meth:`WindowedSweep.apply <orpheus.sn.operators.windowing.WindowedSweep.apply>`
overrides the inherited factor-by-factor body with the substrate's
MOMENT-emit mode (the per-anti-diagonal accumulation of :ref:`Phase 5c
<sn-angular-windowing-in-sweep-accumulation>`), so the
:math:`(N, n_g, n_x, n_y)` full-angular intermediate is **never
materialized** (the ~3× linear peak-memory win). The **inherited**
:meth:`OperatorProduct.apply <orpheus.numerics.operator.OperatorProduct.apply>`
body — ``P.apply(A⁻¹.apply(rhs))``, which *does* materialize the
intermediate — is retained verbatim as the **fuller-view verification
oracle** (the aggressive-retirement oracle exception): the two
evaluations differ only in the ordinate-sum reduction *order*, so they
agree to principled-equivalence within the scale-relative
:math:`4N\varepsilon` bound (measured :math:`1.8\times10^{-16}` on the
product SUT, re-measured from the pre-migration
:math:`2.7\times10^{-16}`). The permanent gate
``test_2d_windowed_product_equals_post_projection`` (renamed/migrated
from ``test_2d_windowed_moments_in_sweep_equal_post_projection`` when the
SUT became the product) is that fused-≡-deforested pin; it is a
*representation*-equivalence pin (procedurally, not structurally,
independent — shared kernel), whose structurally-independent anchors are
the SI ≡ Krylov-full scalar cross-check and the closed-form
:math:`k_\infty` eigenvalue.

.. admonition:: What was tried and rejected — the moment-proxy residual gate
   :class: warning

   An early proposal verified the windowed path with a **moment-proxy
   residual**: compute a residual *in moment space* and read its
   smallness as evidence that ``solve_moments`` inverted its operator. It
   was rejected as **category-confused**. :math:`P` is a coisometry, so
   :math:`P \circ A^{-1}` has no inverse and no round-trip to take a
   residual *of* — a residual gate presumes a solve, but the windowed
   object is a *projection composed with a solve*. Its only honest
   correctness statements are (i) representation-equivalence to the
   deforested oracle and (ii) the structurally-independent scalar anchor.
   Reifying the composition made the confusion **structural**:
   :class:`~orpheus.sn.operators.windowing.WindowedSweep` reports
   ``is_invertible = False`` (its coisometry ``P`` factor is
   non-invertible), so there is no round-trip an inverse-residual gate
   could measure.

At #226 step 2 one transitional wrapper remained: the driver still spoke
``.solve``, so a thin ``_MomentWindowedResolvent`` adapter held the
product and mapped its ``solve`` to the product's ``apply`` (the
``initial_guess`` kwarg accepted for the
:class:`~orpheus.numerics.iteration.SourceIteration` contract and
dropped). **Taxonomy step 3 removed even that** — the SI driver now
consumes inverse-application operators *directly*, so the production SI
holds ``P @ A.inverse()`` with no wrapper. That driver-consumption model
— how the solver builds the inverse and the driver applies it, why an
apply-only step operator is legitimate, and why the adapter could finally
dissolve — is :ref:`inverse-application-driver`.

Cross-references: the reduction :math:`M_{\rm frame}` is
:eq:`harmonic-moment-projection`; the frame's ``analysis`` /
``reconstruction`` faces are the
:class:`~orpheus.numerics.frame.GalerkinFrame`'s
(:ref:`galerkin-projection`); the two-layer "operators, not views" law is
:ref:`operator-algebra`; the reified :math:`M = (L+C-B_{\rm lower})` the
windowed forward may wrap is :ref:`si-gauss-seidel-reification`
(:doc:`discrete_ordinates`); the invertibility-closure and
principled-equivalence framework is ``vv-principles`` § "Bit-identity vs
principled-equivalence".


.. _inverse-application-driver:

The solver builds the inverse; the driver applies it (#226 taxonomy step 3)
===========================================================================

Steps 1–2 reified the inverse *operators* — the schedule-triangular
:class:`~orpheus.sn.operators.sweep_operator.SweepOperator` (step 1, the
WDD sweep as :math:`(L+C)^{-1}`), the reified splitting matrix
:math:`M=(L+C)-B_{\rm lower}` (step 2, :ref:`si-gauss-seidel-reification`
in :doc:`discrete_ordinates`), and the windowed composition
:math:`P\circ A^{-1}` (step 2, :ref:`windowing-retyped` above). **Step 3
retypes the driver–inverse boundary.** The solver (posing) layer builds
the inverse *once*, and the iteration primitive
(:class:`~orpheus.numerics.iteration.SourceIteration`) *applies* it:

.. math::

   \psi_{n+1} \;=\; A_{\rm inv}.\mathrm{apply}
   \Bigl(q_{\rm ext} + \textstyle\sum_i g_i\,\psi_n,\;
   \ \mathrm{initial\_guess}=\psi_n\Bigr),
   \qquad
   A_{\rm inv} \;=\; A^{-1}\ \text{(built by the solver)} .

Concretely :func:`_within_group_si <orpheus.sn.solver._within_group_si>`
does ``step, windowed = _maybe_window(base_resolvent.inverse(), S,
sn_mesh)`` and hands ``step`` to the ``SourceIteration`` constructor as
its **first argument** ``A_inv``. This closes the #226 steps-1–3 arc: the
duck-typed "resolvent" — the object whose ``apply`` and ``solve`` inverted
*different* operators — is **fully dissolved**. Nothing in the driver is a
resolvent any more; there is a forward operator (built by the solver,
invertible by contract) and its inverse (built by the solver, *applied* by
the driver).

This **refines the variadic-driver picture** of
:ref:`bc-extraction-variadic-driver`: the driver's first argument, there
described as "the resolvent it must invert", is — since step 3 — for
:class:`~orpheus.numerics.iteration.SourceIteration` the pre-built inverse
operator it *applies*, while
:class:`~orpheus.numerics.iteration.KrylovAcceleration` keeps the
*forward* :math:`A` for its GMRES matvec (and preconditions with
:math:`A^{-1}`). The homogeneous ``*gains`` bag is unchanged in both.


The seeded-apply family — who threads the start, who ignores it
---------------------------------------------------------------

The inverse family exposes ONE canonical apply signature — the static
contract :class:`~orpheus.numerics.iteration.SupportsSeededApply`:

.. math::

   \mathrm{apply}(\text{rhs},\; *,\; \text{initial\_guess}=\text{None})
   \ \longrightarrow\ \text{codomain} .

The driver threads the previous iterate as ``initial_guess`` on **every**
step, uniformly — the plumbing the curvilinear sweep needs, where the
previous iterate *is* the Morel–Montry Carlson coupled-pole seed (the
sweep reads the level's :math:`\psi(\mu=-1)` from it — the curvilinear
seed obstruction, :ref:`sn-angular-windowing-geometry-restriction`). A
dropped / zeroed / stale seed is a **wrong-fixed-point** bug there, not a
rate change. Members with no use for a start accept the keyword and ignore
it, each documented per type:

.. list-table:: The seeded-apply contract across the inverse family
   :header-rows: 1
   :widths: 34 22 44

   * - Inverse operator
     - ``initial_guess``
     - Why
   * - :class:`~orpheus.sn.operators.sweep_operator.SweepOperator`
       (curvilinear 1-D)
     - **threaded**
     - The M–M Carlson closure seeds the :math:`\mu=-1`
       starting-direction flux from the previous iterate's per-ordinate
       :math:`\psi`; ``apply`` delegates to ``inner.solve(rhs,
       initial_guess=...)`` which reads it.
   * - :class:`~orpheus.sn.operators.sweep_operator.SweepOperator`
       (2-D Cartesian DD)
     - accepted, ignored
     - The diamond-difference wavefront is a direct forward substitution
       down the upwind DAG — no interior-iterate seed. The kwarg threads
       through harmlessly.
   * - :class:`~orpheus.numerics.operator.InverseOperator`
       (value-bearing leaf)
     - accepted, ignored
     - An **exact** pointwise inverse (a division) has no iterative start
       to seed — ``del initial_guess``.
   * - :class:`~orpheus.sn.operators.windowing.WindowedSweep`
       (``P @ A.inverse()``)
     - accepted, ignored
     - The multi-D Cartesian walk has no bulk-seed consumer; a moment
       iterate could not seed the *angular* walk anyway (wrong
       representation). The reflective lag rides the driver's ``B`` /
       ``B_upper`` gain, never this kwarg.

The contract is **static only** (an annotation target), deliberately
**not** ``runtime_checkable`` — mirroring
:class:`~orpheus.numerics.operator.SupportsInverse`. The runtime truth is
the :attr:`~orpheus.numerics.operator.LinearOperator.is_invertible`
property; the Protocol pins the *shape* pyright checks. Threading the seed
**unconditionally** (no per-call decision) is what lets the driver be
shape-agnostic — pinned by the always-on Mode-11 spy
``test_seed_threading_spy.py`` (call *k*'s ``initial_guess`` must equal
call *(k−1)*'s return, by value).


Retiring the signature probe — why it existed, why it could go
--------------------------------------------------------------

Before step 3 the driver could not assume a uniform seed signature,
because it consumed **duck-typed resolvents** whose ``solve`` methods had
heterogeneous signatures (some took ``initial_guess``, some did not). It
therefore probed each one at runtime:

.. code-block:: python

   # RETIRED at #226 step 3
   if _solve_accepts_seed(resolvent):        # inspect.signature(L.solve)
       psi = resolvent.solve(rhs, initial_guess=psi_prev)   # _solve_with_seed
   else:
       psi = resolvent.solve(rhs)

The ``inspect.signature`` probe (``_solve_accepts_seed`` /
``_solve_with_seed``) was a **stringly-typed dispatch on shape** — exactly
the anti-pattern the reification exists to remove. It could retire **only
once the family signature was canonical**: with every inverse advertising
the *same* ``apply(rhs, *, initial_guess=None)``
(:class:`~orpheus.numerics.iteration.SupportsSeededApply`), the driver
threads the keyword unconditionally and the probe has nothing to decide.
Its post-rewire falsifier is the **M-PROBE** mutation (sever the
``initial_guess`` thread at
:meth:`SweepOperator.apply <orpheus.sn.operators.sweep_operator.SweepOperator.apply>`),
which reddens the seed-threading spy in under a second.

Two further gates dissolved with it:

* **The constructor invertibility gate is gone.** The driver's whole
  contract is now a callable ``apply`` — the step operator arrives
  *pre-inverted*, so it never needs a ``solve``. The **invertibility
  obligation moved to the inverse builder**: on a non-invertible leaf
  the inverse simply cannot be constructed (a *structural* leaf declares
  no ``inverse()`` — a static error; a *value-dependent* one raises
  :class:`~orpheus.numerics.operator.NotInvertible`), so "is this a
  faithful inverse of the intended forward?" is discharged where the
  operator is *built*, not re-checked where it is *applied*.
* **An apply-only step operator is legitimate by design.** The windowed
  product ``P @ A.inverse()`` reports ``is_invertible = False`` (its
  :math:`P` factor is a coisometry, :ref:`windowing-retyped`), yet it is
  a valid step operator because a step operator only needs to *apply*.
  Gating on a callable ``apply`` alone — not on invertibility — is what
  makes the apply-only windowed product a first-class step operator
  rather than an illegal state needing an adapter.


The narrowing boundary and the structural resolution (#285)
-----------------------------------------------------------

``SupportsSeededApply`` and ``SupportsInverse`` are not
``runtime_checkable`` (an ``isinstance`` against a Protocol carrying only
a method member is a false-positive machine). The static narrow from "an
invertible operator" to "a seeded-apply operator" therefore lives in
**one** place,
:func:`_seeded_inverse <orpheus.numerics.iteration._seeded_inverse>`:

.. code-block:: python

   def _seeded_inverse(A):                       # numerics/iteration.py
       # runtime precondition: the CALLER's A.is_invertible guard
       return cast(SupportsSeededApply, cast(SupportsInverse, A).inverse())

Both :class:`~orpheus.numerics.iteration.KEigenvalue` (inner
:class:`~orpheus.numerics.iteration.SourceIteration`) and
:class:`~orpheus.numerics.iteration.KrylovAcceleration` (default
preconditioner) route through it, so the cast — and its runtime
precondition, the caller's ``is_invertible`` guard — is single-sourced.

When step 3 shipped, ``_seeded_inverse`` narrowed to the seeded-apply
signature only for the inverses the posing layers here actually reach —
the SN :class:`~orpheus.sn.operators.sweep_operator.SweepOperator` and the
leaf :class:`~orpheus.numerics.operator.InverseOperator`, both of which
carried ``apply(rhs, *, initial_guess=None)`` by per-leaf convention. It
was **not** a claim that *every* ``inverse()`` in the family took the
keyword, and whether that conformance should become **structural** (a
shared mixin, so the signature is *guaranteed* rather than
asserted-per-leaf) or stay per-leaf convention was the open decision
`#285 <https://github.com/deOliveira-R/ORPHEUS/issues/285>`_.

**Step 4 resolved #285 STRUCTURAL.** The wrap-delegate back-half was
extracted into :class:`~orpheus.numerics.operator.InverseWrapMixin`, whose
**abstract** ``apply(x, /, *, initial_guess=None)`` every sibling
inherits — so a new inverse *cannot* forget the keyword: pyright rejects a
kwarg-less override (``reportIncompatibleMethodOverride``,
mutation-verified) and ``ABCMeta`` blocks a sibling that omits ``apply``
entirely. The narrowing ``_seeded_inverse`` performs is therefore now a
cast over an *already-guaranteed* shape rather than a hoped-for one, and
every ``.inverse()`` a posing layer reaches through it carries the seeded
signature by construction.

**Step 5 closed the product residue.** One ``.inverse()`` stayed outside the
family after step 4: a composed
:meth:`OperatorProduct.inverse <orpheus.numerics.operator.OperatorProduct.inverse>`
is a *composition* (:math:`(AB)^{-1}=B^{-1}A^{-1}`), not a wrap-delegate
sibling, and it returned a **raw reversed product** whose positional-only
``apply(x, /)`` raised ``TypeError`` when a driver seeded it. It now returns
:class:`~orpheus.numerics.operator.InverseOperator`\ ``(self)`` — the generic
family member wrapping the product. The action is **bit-identical**
(``apply`` delegates to the product's own ``solve`` = ``b.solve(a.solve(q))``,
exactly the two solves the reversed product composed), but the wrapper adds
the canonical seeded ``apply`` (the #285 ``TypeError`` repro flips to
*accepted* — the solve path never threaded seeds either, so behavior is
unchanged) and **strengthens the involution to object identity**
(``(A@B).inverse().inverse() is (A@B)``, where the reversed-product spelling
rebuilt fresh objects). The factors stay reachable as ``.inner.a`` /
``.inner.b``.

**Two kinds of inverse.** The closure exposes a clean split in what
``.inverse()`` returns across the whole algebra:

.. list-table:: Wrap-delegate vs algebra-closed inverses
   :header-rows: 1
   :widths: 26 34 40

   * - Kind
     - Members
     - What ``.inverse()`` returns
   * - **Wrap-delegate** (the #226 family)
     - :class:`~orpheus.numerics.operator.InverseOperator`,
       :class:`~orpheus.sn.operators.sweep_operator.SweepOperator`,
       :class:`~orpheus.numerics.green_operator.GreenOperator`,
       :class:`~orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator`,
       and now the composed :class:`~orpheus.numerics.operator.OperatorProduct`
     - a thin typed wrapper around the *forward* :math:`A`, carrying the
       canonical seeded ``apply`` — the object *represents* :math:`A^{-1}`
       by inverting :math:`A` on demand.
   * - **Algebra-closed** (first-class forwards)
     - :class:`~orpheus.numerics.operator.PermutationOperator` (inverse *is*
       a permutation), :class:`~orpheus.numerics.operator.IdentityOperator`
       (self-inverse), :class:`~orpheus.numerics.operator.ScaledOperator`
       (:math:`(\alpha L)^{-1}=\alpha^{-1}L^{-1}`)
     - a first-class **forward** operator in the same closed structure — the
       *other* kind of inverse, deliberately left **unwrapped** (a
       permutation's inverse is nothing more exotic than a permutation, so
       wrapping it would only hide its structure).

The wrap-delegate mixin's ``_ForwardT`` bound was also relaxed from
``_InvertibleForward`` to the new minimal ``_WrappedForward`` Protocol
(``domain`` / ``codomain`` / ``apply`` — exactly what the back-half
consumes), which the 4th sibling forced:
:class:`~orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator`
inverts the *materialization* and never touches ``inner.solve`` or
``inner.is_invertible``, exposing the true minimum bound. The mixin and its
enforcement are documented at :ref:`green-operator`; the materialising
sibling, the ``as_matrix`` functor, and the two-error-class boundary at
:ref:`matrix-inverse-operator`.

The two eigenvalue-outer consumers pose the inverse explicitly:

* :class:`~orpheus.numerics.iteration.KEigenvalue` guards
  ``A.is_invertible`` **at construction** (a non-invertible :math:`A`
  raises :class:`~orpheus.numerics.operator.NotInvertible` with a
  domain message, not an ``AttributeError`` mid-iteration) and builds its
  inner :class:`~orpheus.numerics.iteration.SourceIteration` step via
  ``_seeded_inverse(A)``.
* :class:`~orpheus.numerics.iteration.KrylovAcceleration` keeps the
  **forward** :math:`A` (its GMRES matvec is
  :math:`A\cdot - \sum_i g_i\cdot`) and rewired its default-preconditioner
  fallback from the old ``CAP_SOLVE`` probe (with a ``# type: ignore``) to
  the honest ``A.is_invertible`` test + ``_seeded_inverse(A).apply`` — the
  transport-corrected sweep preconditioner (Adams & Larsen 2002 §III).


Verification — the seed spy and the windowed×G-S corner
-------------------------------------------------------

The rewire is pinned by three structural gates (all ``-O``-proof —
``pytest.fail`` / ``np.testing.assert_*``, never a bare ``assert``):

.. list-table:: Step-3 regression gates
   :header-rows: 1
   :widths: 42 58

   * - Gate
     - What it pins
   * - ``test_seed_threading_spy.py``
       (foundation / sentinel; vv Mode 11)
     - Every inner solve's ``initial_guess`` equals the previous iterate's
       return, by value. **Route-invariant across the rewire**: it wraps
       :meth:`InvertibleOperator.solve <orpheus.sn.operators.streaming.InvertibleOperator.solve>`
       — the surface *both* driver generations route through (pre-step-3
       ``resolvent.solve(...)``, post-step-3
       ``SweepOperator.apply`` → ``inner.solve``) — so it was green before
       the rewire and stayed green through it. Teeth: M-SEED-DROP / ZERO /
       STALE (pre-rewire) and **M-PROBE** (post-rewire) all redden it.
   * - ``test_2d_windowed_product_over_gauss_seidel_M_equals_post_projection``
     - The **windowed×G-S corner**: production's 2-D Cartesian default
       ``inner_schedule`` is ``gauss_seidel``, so the driver's *actual*
       step operator is ``P @ M.inverse()`` with
       :math:`M=(L+C)-B_{\rm lower}`. The fused scheduled-walk moment emit
       ≡ the deforested scheduled solve + projection, within scale-relative
       :math:`4N\varepsilon`, with a ``B_lower`` non-degeneracy guard.
       Closes the gap where this corner was pinned only at the
       :math:`\ell=0` / integration level.
   * - ``test_si_single_primitive_contract.py``
     - Both within-group SI paths (eigenvalue inner + fixed-source) build
       the SAME step operator: a
       :class:`~orpheus.sn.operators.sweep_operator.SweepOperator` whose
       ``.inner`` is the
       :class:`~orpheus.sn.operators.streaming.InvertibleOperator` forward
       — the forward identity moved *one level in*, onto ``.inner`` (the
       solver builds the inverse; the driver applies it).

The seed the spy guards is load-bearing **only on curvilinear** meshes: a
simulated seed-drop on the het-2G sphere moves the eigenvalue by
:math:`|\Delta k|\approx3.46\times10^{-2}` (the ``@slow`` value catcher
``test_si_krylov_eigenvalue_equivalence_sphere``, marked
``@catches("ERR-026", "M-SEED-DROP")``). On 2-D Cartesian the direct
wavefront ignores the seed, so the fast net is the *path* spy, not a value
gate. The full tier (SN + numerics + transport) is **2981 passing / 0
real regressions**; the deleted probe also cleared a laundered typed union
from the pyright ratchet (152 → 148).


.. _green-operator:

The Green operator — the preconditioned-splitting sum inverse (#226 step 4)
===========================================================================

Steps 1–3 delivered the *exact/direct* inverse family: the leaf
:class:`~orpheus.numerics.operator.InverseOperator` (a pointwise division)
and the schedule-triangular
:class:`~orpheus.sn.operators.sweep_operator.SweepOperator` (the WDD sweep
as :math:`(L+C)^{-1}`), both invertible in a single pass. **Step 4 adds the
first *iterative* member** — the
:class:`~orpheus.numerics.green_operator.GreenOperator`, returned by
:meth:`OperatorSum.inverse <orpheus.numerics.operator.OperatorSum.inverse>`
— and with it the ``OperatorSum`` invertibility contract that routes to it.
It is the discrete Green's function of a general operator *sum*: the object
whose columns are the flux responses to unit point sources.


The convergent :math:`A`-preconditioned splitting
-------------------------------------------------

A general operator sum has **no operand-wise inverse**: :math:`(A+B)^{-1}`
is not a function of :math:`A^{-1}` and :math:`B^{-1}`
(Sherman–Morrison–Woodbury applies only under low-rank structure), and —
unlike the schedule-triangular ``(L+C)`` — there is no substitution order
that solves it in one pass. What a sum *does* have, when its **leading**
term :math:`A` is invertible, is a convergent **splitting**. Write the sum
as :math:`A - B` (gains carried with their signs) and precondition by
:math:`A^{-1}`:

.. math::
   :label: green-neumann-series

   (A - B)^{-1}\,q
   \;=\; \bigl(I - A^{-1}B\bigr)^{-1} A^{-1}\,q
   \;=\; \sum_{k=0}^{\infty} \bigl(A^{-1}B\bigr)^{k}\,A^{-1}\,q .

.. vv-status: green-neumann-series documented

The right-hand side is the **Neumann series** of the :math:`A`-preconditioned
splitting. In transport it *is* the **multiple-scattering expansion**:
:math:`A^{-1}q` is the uncollided flux, :math:`(A^{-1}S)\,A^{-1}q` the
once-rescattered contribution, :math:`(A^{-1}S)^{k}A^{-1}q` the
:math:`k`-times-rescattered contribution ([LewisMiller1984]_ §3.2;
[AdamsLarsen2002]_ §II). Its partial sums are **exactly** the Richardson /
source-iteration iterates started from zero,

.. math::
   :label: green-splitting-iteration

   x_{n+1} \;=\; A^{-1}\bigl(q + B\,x_n\bigr),
   \qquad x_0 = 0
   \;\Longrightarrow\;
   x_{n} = \sum_{k=0}^{n-1}\bigl(A^{-1}B\bigr)^{k}A^{-1}q ,

.. vv-status: green-splitting-iteration documented

which converge iff the iteration matrix is a contraction,
:math:`\rho(A^{-1}B) < 1`. For the within-group transport loss
:math:`A_{\rm loss} = (L+C) - S` this is the physical **scattering-ratio**
bound

.. math::

   \rho\bigl((L+C)^{-1}S\bigr)
   \;\le\; \max_{\rm cell}\ \frac{\Sigma_s}{\Sigma_t} \;=\; c \;<\; 1 ,

guaranteed below unity for any absorbing (physical) medium
([AdamsLarsen2002]_ §II: :math:`\rho = c`). The convergence is thus a
material property, not a numerical accident — and it is precisely why the
sum has an inverse *operator* at all.


The name is earned — G-Neumann, G-reciprocity, G-kernel
-------------------------------------------------------

A subclass name in this family is a **promise backed by a test** (taxonomy
§13): a distinguishing invariant a *bare*
:class:`~orpheus.numerics.operator.InverseOperator` does not automatically
have. Round-trip alone (:math:`A^{-1}(Ax)=x`) earns only the generic name;
``GreenOperator`` must be *tested* to carry the Green's-function structure.

.. list-table:: The three name-earning invariants of ``GreenOperator``
   :header-rows: 1
   :widths: 22 78

   * - Invariant
     - What it asserts (and why a generic :math:`A^{-1}` fails it)
   * - **G-Neumann**
     - The partial sums of :eq:`green-neumann-series` equal
       ``green.apply(q)`` (the multiple-scattering expansion). A generic
       :math:`A^{-1}` has **no splitting** to satisfy — this is the
       *distinguishing* invariant. Pinned with the exact geometric decay of
       the split (see :ref:`green-verification`).
   * - **G-reciprocity**
     - The Green's-function reciprocity theorem in the **Euclidean** inner
       product, :math:`\langle\phi_2, G\phi_1\rangle = \langle
       G^{\mathsf T}\phi_2, \phi_1\rangle`, where :math:`G^{\mathsf T}` is
       the Green built over the **transposed operands** (:math:`A^{\mathsf T}
       = A^{\mathsf T}_{\rm lead} - B^{\mathsf T}`). It is the *cheap* proof
       (no second dense oracle) that the split derivation is correct for a
       *different* operand configuration.
   * - **G-kernel**
     - ``green.apply(δ_j)`` is column :math:`j` of :math:`(A-B)^{-1}` — the
       flux response to a unit point source at :math:`j`. Folded into the
       anchor's input set (the :math:`\delta_j` basis) rather than gated
       separately, since it *is* the forward anchor evaluated on unit
       vectors.

.. important::

   **G-reciprocity uses the Euclidean transpose, not the metric adjoint.**
   :math:`G^{\mathsf T}` here is the plain transpose under the :math:`L^2`
   pairing, built manually from the transposed operands — **not** the
   ``.H`` Hilbert adjoint :math:`G^{\dagger} = \mathcal G^{-1}
   G^{\mathsf T}\mathcal G` carrying the angular Gram metric. The
   adjoint-inverse axis (``G.H == (A-B).H.inverse()`` — free at the object
   level on the iterative branch, but blocked on the SN preconditioner's
   as-yet-unbuilt multi-D transpose sweep) is the separate
   `#280 <https://github.com/deOliveira-R/ORPHEUS/issues/280>`_ family and
   is deliberately *not* promised here.

The family is keyed by **which mathematical object the inverse is, not by
the algorithm that realizes it**. A Richardson-realized Green and a
GMRES-realized Green are the *same* Green operator — which is why
``KrylovInverseOperator`` was **rejected** as a sibling name: "Krylov"
names the orthogonal *realization* axis, not the object.


Green wraps the driver — it re-implements nothing
-------------------------------------------------

``GreenOperator`` is a thin wrapper over a driver, not a new solver
(taxonomy §11.2, Pattern 5). Its application engine is
:class:`~orpheus.numerics.iteration.SourceIteration` — the **same**
Richardson driver the SN solver and
:class:`~orpheus.numerics.iteration.KEigenvalue` consume directly. At
construction it derives the splitting **once, from the sum's structure**:

* the **left-spine head** becomes the preconditioner :math:`A^{-1}`,
  through *its own* structure-keyed ``.inverse()`` — the WDD
  :class:`~orpheus.sn.operators.sweep_operator.SweepOperator` for an
  ``(L+C)`` head, a leaf
  :class:`~orpheus.numerics.operator.InverseOperator` for a value-bearing
  head;
* every remaining term becomes a **negated gain** of the driver (the
  ``A - S`` spelling arrives as ``ScaledOperator(-1, S)``; the gain is
  :math:`S` itself, un-wrapped so the driver holds the *named* operator).

It then hands these to ``SourceIteration`` and re-implements no iteration
math of its own. The left spine is flattened by walking **exact**
``OperatorSum`` nodes only: the fused
:class:`~orpheus.sn.operators.streaming.InvertibleOperator` is a structural
*leaf* (its sum-ness is an MRO fact, its identity a fused operator with a
direct inverse), so ``((L+C) - S)`` flattens to preconditioner ``(L+C)`` +
gain ``[S]`` — never dissolving the ``(L+C)`` into its own summands. A
GMRES-realized Green (with
:class:`~orpheus.numerics.iteration.KrylovAcceleration` as the engine) is a
future realization *strategy* of this same object, not a sibling type.


The ordering ruling — four edges of one canonical order
-------------------------------------------------------

Operand **spelling** selects the algorithm — the #261 canonical-ordering
rule, extended to the sum inverse. The four edges:

.. list-table:: How a sum's spelling routes its inverse
   :header-rows: 1
   :widths: 20 26 54

   * - Spelling
     - ``.inverse()`` builds
     - Behaviour
   * - ``L + C``
     - :class:`~orpheus.sn.operators.sweep_operator.SweepOperator`
     - The fusion dispatch on
       :meth:`StreamingOperator.__add__ <orpheus.sn.operators.streaming.StreamingOperator.__add__>`
       (streaming.py:510, "the canonical and only ordering") returns the
       :class:`~orpheus.sn.operators.streaming.InvertibleOperator`
       specialisation, whose ``.inverse()`` override (→ the direct sweep)
       **shadows the generic Green by MRO** (type-as-structure, §11.1).
   * - ``A_loss = (L+C) - S``
     - :class:`~orpheus.numerics.green_operator.GreenOperator`
     - Leading ``(L+C)`` invertible → the physical splitting: preconditioner
       the sweep, gain :math:`S`. Converges (:math:`c<1`); does **not**
       raise.
   * - ``C + L``
     - :class:`~orpheus.numerics.green_operator.GreenOperator`
     - A *legal* spelling whose leading term ``C`` happens to be invertible,
       so a Green **constructs** — the algebra cannot read
       :math:`\rho(C^{-1}L) > 1`. Its collision-preconditioned Richardson
       **diverges**, and ``apply`` raises
       :class:`~orpheus.numerics.green_operator.ConvergenceFailure`
       **loudly**. Same math as ``L + C``, different algorithm by spelling —
       never a silent wrong answer (Cardinal Rule 1).
   * - ``(-S) + A``
     - *(refused)*
     - The left-spine head :math:`-S` is not invertible; the factory raises
       :class:`~orpheus.numerics.operator.NotInvertible` at
       **construction**, naming the canonical ordering (spell the invertible
       operator first, ``A - S``).

The keying predicate is honest about its scope.
:attr:`OperatorSum.is_invertible <orpheus.numerics.operator.OperatorSum.is_invertible>`
returns ``self.a.is_invertible`` — the *leading* term — and therefore reads
"**leading-term-preconditionable at this operand order**", **not**
spelling-independent mathematical invertibility (verification spec §18.B).
For ``(-S) + A`` the operator :math:`A - S` *is* mathematically invertible,
yet the predicate reports ``False`` because :math:`-S` is spelled first.
This is acceptable precisely because the #261 rule already makes operand
order semantically load-bearing, no production consumer relies on the
spelling-independent meaning, and the refusal is **loud**. Since carve
P4 there is no parallel ``CAP_SOLVE`` tag to keep in lockstep:
``is_invertible`` is the *single* advertisement, and the keystone v2
faithfulness contract (:ref:`capability-set-semantics`) pins that
``.inverse()`` returns exactly when ``is_invertible`` is ``True``.


The promise — the TRUE residual, driven not merely checked
----------------------------------------------------------

``GreenOperator`` promises the **converged** :math:`A^{-1}q` or a loud
:class:`~orpheus.numerics.green_operator.ConvergenceFailure` — never a
silent partial iterate. The subtlety is *which residual* defines
"converged". :class:`~orpheus.numerics.iteration.SourceIteration` stops on
the iterate **increment** :math:`\lVert\Delta\psi\rVert / \lVert\psi\rVert
< {\rm tol}`, which **understates** the true equation error by the factor
:math:`\rho/(1-\rho)` (numerical-bug-signatures Signature 9, ρ-blind
stopping). As :math:`\rho \to 1` an increment-converged iterate can sit
orders of magnitude off the equation. The promise is therefore read on the
**true relative residual**:

.. math::
   :label: green-true-residual

   \frac{\bigl\lVert (A - B)\,\psi - q \bigr\rVert}{\lVert q \rVert}
   \;<\; {\rm tol} .

.. vv-status: green-true-residual documented

Crucially the promise is **driven, not merely checked**. A check-only design
(raise whenever :eq:`green-true-residual` fails at increment-stop) would
falsely raise for **every** split with :math:`\rho > 1/2` — i.e. most
physical scattering ratios — because increment-stop delivers only
:math:`\rho/(1-\rho)\cdot{\rm tol}` there.
:meth:`GreenOperator.apply <orpheus.numerics.green_operator.GreenOperator.apply>`
instead runs a **refinement loop**: after each increment-stopped driver
call it measures :eq:`green-true-residual` and, if unmet with budget
remaining, re-seeds the driver with its own iterate (steps accumulate
against one total ``max_iter``). The driver stays the sole iteration engine;
the loop is tolerance bookkeeping only, at one extra forward matvec per
check.

The refinement's terminal ``raise`` also has to be **NaN-safe**, because a
hard-divergent split (the ``C + L`` trap) produces two distinct
floating-point failure shapes — both found by the divergence tooth
(2026-07-02):

#. the iterate itself overflows, so the increment is ``nan``; and
#. — sharper — the driver's stopping **denominator** :math:`\lVert\psi\rVert`
   overflows to ``inf`` one step *before* the numerator, so the driver
   "converges" at ``increment = finite/inf = 0.0`` onto a ~\ :math:`10^{154}`
   garbage iterate.

Both are caught because the promise test reads the true residual of the
returned iterate (huge or non-finite in either shape), never the driver's
increment. This is exactly the "what was tried and failed" material a naive
increment-check would have shipped as a silent wrong answer.


Green versus the k-eigenvalue inner — normal form versus inexact relaxation
---------------------------------------------------------------------------

``K = A_loss.inverse() @ F`` is the **normal form** of the k-eigenvalue
problem: the eigen-operator whose dominant eigenvalue is :math:`k`. But
production power iteration **deliberately does not** build this exact
``GreenOperator`` — it keeps consuming
:class:`~orpheus.numerics.iteration.SourceIteration` directly, as a
**warm-started, budget-bounded inner relaxation**. That partial convergence
is *by design* in nested iteration (classic inexact power iteration): the
inner need only be converged enough that the *outer* dominant-eigenvector
direction is accurate, and warm-starting from the previous outer iterate
amortises the cost. This is a **different contract** from Green's
converged-or-raise promise, and conflating them would either over-solve the
inner (wasted matvecs) or import Green's hard raise into a loop that
legitimately runs partial solves. So
:class:`~orpheus.numerics.iteration.KEigenvalue` was examined and
**deliberately not rewired** to ``A_loss.inverse() @ F``; the normal form is
the *concept*, the inexact inner is its production *realization*.

Green's consumers today are the invariant gates (see
:ref:`green-verification`). Production consumers arrive later: the #200
preconditioner algebra, a diffusion ``.inverse()`` realized as a
CG-preconditioned Green (the taxonomy's *negative control* — an iterative
inverse with no sweep), and future explicit normal-form spellings of the
k-problem.


The wrap-delegate mixin — extracted at the third sibling
--------------------------------------------------------

Every inverse in this family is a thin typed wrapper around its own forward
operator :math:`A`: ``apply`` realizes :math:`A^{-1}` by the sibling's
algorithm, and *everything else is delegation*. That back-half —
``is_invertible = True`` with ``solve`` = the forward matvec
``inner.apply`` (solving :math:`A^{-1}y = b` *is* applying :math:`A`),
the domain↔codomain swap
(an inverse maps the forward's codomain back to its domain), and the
object-identity involution ``inverse() → inner``
(:math:`(A^{-1})^{-1} = A`) — was carried byte-identically by
:class:`~orpheus.numerics.operator.InverseOperator` and
:class:`~orpheus.sn.operators.sweep_operator.SweepOperator` as documented
twins. ``GreenOperator`` is the **third** sibling, and it fired the
extraction trigger both twins recorded (defer-until-≥2, extract at 3): the
shared mechanism lifted into
:class:`~orpheus.numerics.operator.InverseWrapMixin`. The forward-operand
Protocol was renamed ``_SolveBackedLeaf`` → ``_InvertibleForward`` to match
its now-general role as the mixin's ``_ForwardT`` bound. Each sibling keeps
exactly three things of its own: the constructor **guard** (what makes its
``inner`` invertible — a value check, a type, a derivable splitting), the
``apply`` **body** (the inversion algorithm), and ``__repr__``.

The mixin also carries the **abstract** seeded-apply signature ``apply(x,
/, *, initial_guess=None)`` — the structural resolution of #285 discussed at
:ref:`inverse-application-driver`: a new sibling *cannot* forget the keyword
(pyright rejects a kwarg-less override; ``ABCMeta`` blocks a missing one).
The adjoint axis is **not** part of the back-half — ``is_adjointable`` /
``.H`` stay at the base defaults, deferred to the #280 family.


.. _green-verification:

Verification
------------

``GreenOperator`` is the **first** inverse in the family with **no legacy
``.solve`` to inherit from** — the sum was not invertible before step 4, so
there is no bit-identity twin to ride. Its correctness therefore rests
**entirely on structurally-independent anchors**, never on inheritance.
Every Green value gate is a foundation / flux-shape claim against such an
anchor: **no eigenvalue claim, no MMS reference** (an iterative sum inverse
is not an eigenvalue solver, and MMS is source-driven — neither pillar
applies here). The ``inverse().apply ≡ solve`` equivalence that anchored
steps 1–3 is a **tautology** for the sum (since carve P4 the generic sum
carries **no** ``solve`` at all — ``inverse().apply`` *is* its only
inverse action, so no second spelling exists to compare) and is
deliberately excluded as evidence.

.. list-table:: Step-4 verification gates (all ``@pytest.mark.foundation``)
   :header-rows: 1
   :widths: 42 58

   * - Gate
     - What it pins
   * - ``test_green_operator.py`` (L0, dense split ``A = D − αP``)
     - **G-I1** round-trip both ways at driver tol + the **dense-LU anchor**
       (:func:`numpy.linalg.solve` of the materialized sum), including the
       :math:`\delta_j` basis (the G-kernel fold). **G-Neumann** partial
       sums converge to ``green.apply(q)`` with the **exact** 4-cycle decay
       ratio :math:`\rho = \alpha/(\prod d_i)^{1/4}` (the permutation's
       spectrum makes the pin exact, not a fuzzy band). **G-reciprocity**
       via the transposed-operand Green (non-symmetric ``A``, so not
       vacuous). Divergence + near-critical (:math:`\rho=0.99`) raises, each
       with a convergent control.
   * - ``test_green_operator_sn.py`` (L1, het 2G **vacuum** slab)
     - **G-Neumann-L1** on the real operators: :math:`\sum_k ((L+C)^{-1}S)^k
       (L+C)^{-1}q` → ``green.apply(q)`` with the geometric tail of the
       physical scattering ratio. The anchor is a **trace-consistent
       manufactured** pair (``x_tc = (L+C)⁻¹(random)``, ``q = A_loss·x_tc``),
       which resolves the #284 source-subspace caveat — sweep inverses are
       exact on the source subspace, so the exact solution *is* ``x_tc`` with
       no dense-LU trace mismatch. Plus **driver bit-identity** vs the
       hand-built ``SourceIteration(sweep, S)``, and the four ordering-ruling
       edges end-to-end.

The config discipline is Mode-9: the L1 slab is **heterogeneous, ≥2-group,
vacuum** — a reflective isotropic box would null the
streaming↔scattering redistribution the Neumann series expands, and 1-group
is blind to a scattering-matrix transpose (Mode 6). The teeth are
**14 mutations verified** (2026-07-02): 12 bite their named gates
(sign/swap/flatten/tol/increment/seed/order/…), and 2 are designed-green
controls — the always-wrap ``_negated`` no-op (proving the gain-unwrap is
pure deforestation) and the **M-GRN-SEED blindness proof** (a dropped
driver-start reddens Green's own Mode-11 seed spy while the landed step-3
spy and every value gate stay green, proving the step-3 spy is structurally
blind to Green's driver-start threading). The pyright ratchet held exactly
at baseline; the Sphinx build stayed ``-W`` clean.


.. _matrix-inverse-operator:

The materialising functor and the dense direct inverse (#226 step 5)
====================================================================

Steps 1–4 built the inverse *operators* — the exact leaf
:class:`~orpheus.numerics.operator.InverseOperator`, the
schedule-triangular
:class:`~orpheus.sn.operators.sweep_operator.SweepOperator`, and the
iterative :class:`~orpheus.numerics.green_operator.GreenOperator` — each a
matrix-free wrapper that realizes :math:`A^{-1}` by *applying* it. **Step 5
adds the** *materialising* **family**: the serialization boundary
:meth:`~orpheus.numerics.operator.LinearOperator.as_matrix` promoted to a
universal ``LinearOperator`` base method, and the fourth
:class:`~orpheus.numerics.operator.InverseWrapMixin` sibling
:class:`~orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator` —
the inverse of a **structureless small** operator, obtained by
materializing :math:`[A]`, factoring it once (LU), and back-solving every
:meth:`apply`. The step also **closes #285** for composed inverses (the product residue
retired in the #285 structural-resolution section above) and retires the
homogeneous solver's ``_as_dense`` prototype into the new base method.


``as_matrix`` — the functor out of the operator category
--------------------------------------------------------

The inverse, adjoint, and composition maps of :ref:`operator-algebra` are
all **endofunctors** :math:`\mathrm{Op}\to\mathrm{Op}`: each takes an
operator to another operator. ``as_matrix`` is the **fourth arrow** of the
taxonomy (§2) — the functor *out* of the operator category,
:math:`\mathrm{Op}\to\mathrm{Mat}`, the serialization boundary that leaves
the matrix-free world entirely. It is realized as the **apply-to-basis**
loop: column :math:`j` is the operator applied to the :math:`j`-th basis
element,

.. math::
   :label: matrix-functor-out

   [A]_{:,j} \;=\; \operatorname{ravel}\bigl(A\,e_j\bigr),
   \qquad
   e_j \;=\; \operatorname{unravel}(\delta_j),
   \qquad
   j = 0,\dots,\textstyle\prod_k b_k - 1,

.. vv-status: matrix-functor-out documented

with basis elements enumerated in **C-order** over the carrier shape
:math:`(b_0,\dots,b_{d-1})` and each output raveled the same way, so
:math:`[A]\,x_{\rm flat} = (A\,x)_{\rm flat}` exactly and — for a
group-leading :math:`(n_g,1)` carrier — column :math:`j` is the response to
a unit source in group :math:`j`. The matrix is
:math:`(\prod\text{out shape})\times(\prod_k b_k)`: **rectangular whenever
the operator is not endomorphic**, because the output dimension emerges
from :meth:`apply` itself, never from declared metadata (a coisometry
materializes as a genuinely wider-than-tall block, :ref:`windowing-retyped`).

This default body is the **promoted** ``_as_dense`` apply-to-basis loop
that lived privately in ``homogeneous/solver.py`` until step 5. Its
promotion to the base class is the Cardinal-Rule-2 move: the one place that
turned an operator into a dense matrix becomes the method *every* operator
inherits, and the homogeneous solver drops its private copy
(:ref:`matrix-inverse-consumers`).


The resolution rule and the two error classes
---------------------------------------------

Two questions must be answered before the loop runs: *what shape* are the
basis elements, and *is the materialization affordable*. They raise **two
deliberately distinct exception classes** (spec §27.C), because they are
different kinds of failure.

**Basis-shape resolution** is single-sourced in ``_resolve_basis_shape``
(shared with the eager
:class:`~orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator`
constructor, which must know the resolved shape to reshape solutions back
into carriers — one source, no reciprocal drift). The rule has three arms:

#. an explicit ``basis_shape`` wins;
#. else the operator's own
   :attr:`~orpheus.numerics.operator.LinearOperator.domain` supplies
   ``domain.shape``;
#. else — **no domain and no explicit shape** — the request is *ill-posed*:
   the basis cannot be derived at all, and the method raises ``ValueError``
   naming both remedies (construct the operator with a space, or pass an
   explicit ``basis_shape=(n_g, 1)``).

**The size gate** is orthogonal. When the resolved dimension
:math:`n=\prod_k b_k` exceeds ``max_dimension``, the method raises
:class:`~orpheus.numerics.operator.MatrixTooLarge` — a *well-posed* request
that this environment *refuses on resource grounds*.

The two must never be caught with one loose ``except (ValueError,
RuntimeError)``: a bug that collapsed the two would pass such a test. The
boundary gates pin each separately, and their ``pytest.raises`` matches are
class-discriminating — a Pattern-4 illegal-states check, "un-materializable
*as posed*" :math:`\neq` "too big to materialize *here*".


The size gate and the dense-vs-sparse ruling
--------------------------------------------

:class:`~orpheus.numerics.operator.MatrixTooLarge` is a **RuntimeError**,
not a ``ValueError`` — and there is deliberately **no** ``is_materializable``
predicate beside
:attr:`~orpheus.numerics.operator.LinearOperator.is_invertible`. The reason
is the taxonomy §17 A2 distinction: ``as_matrix`` is a **total functor**.
*Every* finite-dimensional linear operator *has* a matrix; nothing about
its structure or values can make it "un-materializable". What the gate
guards is a **resource effect** — materializing commits :math:`O(n^2)`
memory and :math:`n` applications — so the refusal is a runtime resource
precheck (it reads nothing but a dimension), *not* a structural restriction
(``is_invertible`` / ``is_adjointable`` read the operator's structure *and*
values). The default budget is ``max_dimension = 4096`` — a :math:`4096^2`
float64 is 134 MB and 4096 applications, generous for every
dense-by-construction consumer (0-D energy spectra, CP :math:`[P]`) and
prohibitive for a meshed SN full-field operator by design. It is a per-call
knob: ``try: A.as_matrix() except MatrixTooLarge: <iterative path>``, or
raise the budget for one call.

**The return is a dense** :class:`numpy.ndarray`. The dense-vs-sparse
return keying was the one open thread (W4) that taxonomy §11.4 explicitly
**deferred to step 5**, and the ruling is: *keyed by the operator's
structural override, with dense the only realization built until a sparse
consumer exists*. A structured operator MAY override ``as_matrix`` with a
direct assembly — the future per-octant sparse-triangular streaming
assembly noted at ``sweep_graph.py:66`` is exactly such an override, and it
is **deferred with its 3-D consumer** (there is no sparse consumer today,
so building the sparse path now would be a primitive with no product —
defer-until-consumer, Cardinal Rule 2).
:class:`~orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator` is
the one override that ships in step 5, and it collapses to a single batched
LU back-solve (below), not a sparse form.

.. note:: **Landed (stencil-assembly 2b, 2026-07-04).** The deferred
   sparse-assembly override arrived. Not yet the *octant-batched 3-D*
   form — that vectorization is still deferred with its 3-D consumer —
   but the assembly axis itself: structured leaves now emit a
   :class:`~orpheus.numerics.assembled_operator.SparseAssembledOperator`,
   and ``as_matrix`` **delegates** to the densified emission when the
   operator is :func:`~orpheus.numerics.operator.assemblable` (the
   probing loop retained as the fallback and fuller-view oracle). The
   dense-vs-sparse keying above is now *realised* rather than deferred —
   the key is exactly the ``is_assemblable`` predicate. See
   :ref:`operator-algebra-assembly-axis`.


``MatrixInverseOperator`` — the dense direct inverse
----------------------------------------------------

``A.inverse()`` is a **structure-keyed factory** (taxonomy §13): the
concrete subclass names the mathematical *object* the inverse is, never the
algorithm. A schedule-triangular forward returns the direct-substitution
:class:`~orpheus.sn.operators.sweep_operator.SweepOperator`; a general sum
with an invertible leading term returns the preconditioned-splitting
:class:`~orpheus.numerics.green_operator.GreenOperator`; a value-bearing
leaf returns the pointwise
:class:`~orpheus.numerics.operator.InverseOperator`.
:class:`~orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator` is
the inverse of a **structureless small** operator — the 0-D energy
spectrum, a CP collision-probability block, any composition whose only
exploitable property is that it *fits*. Its algorithm is the honest
consequence of "no structure to exploit, but small enough to hold":
materialize :math:`[A]` once, ``lu_factor`` it once, and every
:meth:`apply` is a direct ``lu_solve`` back-solve; its
:meth:`~orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator.as_matrix`
override collapses the base loop to one **batched** ``lu_solve(lu, I)``
against the same factors.

It is the fourth
:class:`~orpheus.numerics.operator.InverseWrapMixin` sibling, so its entire
back-half — ``is_invertible = True``, the domain↔codomain swap,
``solve`` = the forward matvec ``inner.apply``, and the object-identity
involution ``inverse() → inner`` — is inherited. It keeps only the family's
three per-sibling pieces: its constructor guard, its ``apply`` body, and
``__repr__``. Its home is a leaf module (the ``green_operator.py``
placement precedent), and — *unlike* Green — nothing in ``operator.py``
routes back to it: **no automatic** ``.inverse()`` **returns this type**,
so there is no late-import seam at all, and (step 5 scope) it is **direct
construction only**, never a factory-dispatch target. The
``A.inverse() → MatrixInverseOperator`` routing and the normal-form
spelling ``K = MatrixInverseOperator(loss) @ F`` land later (task #138 / CP).


The name is earned — M-materialise and M-direct at the precision grain
----------------------------------------------------------------------

A subclass name in this family is a **promise backed by a test**
(taxonomy §13): a distinguishing invariant a bare
:class:`~orpheus.numerics.operator.InverseOperator` does not automatically
have. For ``MatrixInverseOperator`` the promise is the **Matrix** invariant
in two faces — M-materialise (the explicit inverse, both ways) and M-direct
(the true residual):

.. math::
   :label: matrix-inverse-materialise

   \bigl\lVert\,[A^{-1}]\,[A] - I\,\bigr\rVert
   \;=\;
   \bigl\lVert\,[A]\,[A^{-1}] - I\,\bigr\rVert
   \;\le\; K\,\varepsilon_{\rm mach}\,\kappa(A)
   \qquad\text{(M-materialise, both ways)} ,

.. vv-status: matrix-inverse-materialise documented

.. math::
   :label: matrix-inverse-direct-residual

   \frac{\bigl\lVert A\,(A^{-1}q) - q\bigr\rVert}{\lVert q\rVert}
   \;\le\; K\,\varepsilon_{\rm mach}\,\kappa(A)
   \qquad\text{(M-direct, seed-independent)} .

.. vv-status: matrix-inverse-direct-residual documented

**What earns the name is the precision GRAIN, not the existence of a
matrix** — and this **supersedes the taxonomy §13 M-row parenthetical**
("a matrix-free inverse *cannot* satisfy M-materialise, having no
``as_matrix``"). That parenthetical is now **false**: step 5 made
``as_matrix`` universal, so an *iterative*
:class:`~orpheus.numerics.green_operator.GreenOperator` **also** satisfies
:math:`[\,\text{green}\,]\,[A]\approx I` — each column ``green.apply(e_j)``
:math:`\approx A^{-1}e_j`. "No ``as_matrix``" no longer distinguishes the
direct inverse; the **grain of the approximation** does:

.. list-table:: Why the grain is the name-earner (spec §27.A)
   :header-rows: 1
   :widths: 24 38 38

   * - Face
     - ``MatrixInverseOperator``
     - an iterative Green over the same :math:`A`
   * - :math:`[A^{-1}]\,[A]-I`
     - one batched ``lu_solve`` on the identity against the SAME
       factorization — **machine·cond**, no second realization, no
       iteration floor
     - :math:`n` iterative solves, each to **driver tolerance**
       (:math:`\sim 10^{-8}`) — floors at its stopping tol
   * - :math:`\lVert A(A^{-1}q)-q\rVert`
     - a direct solve — **machine·cond** true residual
     - floors at the driver tolerance, never machine grain

So M-materialise and M-direct are gated at **machine·cond** grain
(:math:`\text{atol}=K\,\varepsilon_{\rm mach}\,\kappa`), and the gate carries
an explicit **contrast**: the *same* :math:`A`, wrapped as a
:class:`~orpheus.numerics.green_operator.GreenOperator`, meets only
driver-tol — proving the invariant *distinguishing*, not merely satisfied
(the §21 discipline "a bare invariant a generic inverse also satisfies is
not name-earning"). M-direct is also the **seed-independence** face: an
exact direct inverse has nothing to seed, so the canonical
``initial_guess`` keyword is accepted and *ignored*, and the result is
bit-identical under any seed — contrast the sweep's Carlson closure and
Green's splitting start, which *consume* it.


Values, not structure — the guard difference and the witness
------------------------------------------------------------

The constructor deliberately **does not consult** ``inner.is_invertible``.
That predicate is *structural* — it reads the operand tree, and for a sum it
reports "leading-term-preconditionable at this spelling"
(:ref:`green-operator`), a property of the *spelling*, not of the matrix.
``MatrixInverseOperator`` reads **values**: it materializes :math:`[A]` and
factors it, and a matrix is either numerically invertible or it is not,
regardless of how the operator that produced it was spelled. The guards
that remain are the honest value-level ones, all raised at **construction**
(this module family's composition-time-not-call-time principle), never at
apply time:

.. list-table:: The three construction-time guards
   :header-rows: 1
   :widths: 22 78

   * - Guard
     - Behaviour
   * - **Size**
     - :class:`~orpheus.numerics.operator.MatrixTooLarge` propagates from
       the eager materialization (the size gate above).
   * - **Squareness**
     - a rectangular materialization has no two-sided inverse — a
       ``ValueError`` in domain language ("M-materialise is unsatisfiable").
   * - **Exact singularity**
     - a zero LU pivot raises :class:`numpy.linalg.LinAlgError`. scipy's
       *own* singularity signal is only a ``LinAlgWarning`` (``getrf`` info
       :math:`>0`), which would let a zero pivot flow into ``inf`` / ``nan``
       back-solves — the constructor **silences that warning and raises the
       loud error** from the U-diagonal instead (Cardinal Rule 1: fail at
       construction, never return a non-inverse). **Near**-singularity is
       *not* refused — it is priced into the M-direct :math:`\kappa(A)`
       bound.

**The witness.** The values-vs-structure difference is provable, and it is
the taxonomy §3 ``strategy=`` override seam realized *honestly* — not a flag
on ``.inverse()`` but an explicit construction by a consumer who knows the
problem is small (**the type IS the strategy choice**). The witness is a sum
whose **leading term is not invertible**: the canonical motivating form is
the SN :math:`(-S)+(L+C)`, which
:class:`~orpheus.numerics.green_operator.GreenOperator` **refuses at
construction** (left-spine head :math:`-S` non-invertible, canonical-ordering
:class:`~orpheus.numerics.operator.NotInvertible`) yet which materializes
to a perfectly invertible matrix. Because the real :math:`(-S)+(L+C)` is a
``FullField`` carrier — **out of ``as_matrix``'s ndarray scope** (the
honest-scope note below) — the gate proves the identical refusal/inversion
asymmetry on the realizable ndarray analog :math:`(-S_{\rm ao})+D` (an
apply-only leaf :math:`S_{\rm ao}`, ``is_invertible = False``, plus a
diagonal :math:`D`): ``.inverse()`` refuses it, while
``MatrixInverseOperator`` materializes :math:`D-S_{\rm ao}` and inverts it,
anchored against :math:`\mathrm{np.linalg.solve}(D-S_{\rm ao},\,q)`. The
*structure* — a leading-non-invertible sum that is nonetheless an invertible
matrix — is what the witness proves; the physics is incidental.

.. note::

   **Honest scope — ``as_matrix`` serves ndarray carriers only.** The
   apply-to-basis loop builds bare-ndarray :math:`e_j`; a **typed-carrier**
   (``FullField`` SN composite) operator's ``apply`` needs a ``FullField``
   basis vector and sits far above any sane size gate, so it stays
   matrix-free. Step 5 gates ndarray-carrier operators only (diagonal /
   permutation / the model-shared energy leaves / small compositions /
   dense test operators); the ``FullField`` :math:`(L+C)` and
   :math:`(-S)+(L+C)` materialization is a future carve (its 3-D
   sparse-triangular sibling is the deferred ``sweep_graph.py:66``
   override). This is exactly why the witness above uses the ndarray analog.


.. _matrix-inverse-consumers:

The consumer ruling — the latent normal form
--------------------------------------------

Like :class:`~orpheus.numerics.green_operator.GreenOperator` at step 4,
``MatrixInverseOperator`` **shipped** at step 5 verified but not yet wired as
a production spelling — its value gates constructed it directly, and the
factory routing waited. Taxonomy **step 5b closed that loop**: the homogeneous
solver is now the **first production consumer** (below). What still waits is
the *factory routing* — an automatic ``.inverse()`` that returns this type —
which homogeneous deliberately **bypasses** (see the ``direct_eigenvalue``
bullet). The consumers, from latent to landed:

* **The retired ``_as_dense``.** The homogeneous solver's private
  apply-to-basis loop was retired into
  :meth:`~orpheus.numerics.operator.LinearOperator.as_matrix` at step 5, so
  the materialization of the loss operator :math:`\mathbf A = C - K_{\rm iso}`
  and the fission dyad :math:`\mathbf F = \chi\otimes\nu\Sigma_f` both flow
  through the one promoted base method rather than a bespoke loop. The output
  stayed **byte-identical** through that retirement (same basis columns, same
  C-order, same eigen call), and the landed SymPy :math:`k_\infty` pins stayed
  green untouched. (Step 5b then re-spelled the whole eigensolve as
  ``K = MatrixInverseOperator(loss) @ production``: the loss materialization
  now happens *inside* the inverse operator's constructor, and the product's
  own ``as_matrix`` forms :math:`[\mathbf K] = \mathbf A^{-1}\mathbf F`.) See
  :ref:`theory-homogeneous`.
* **``direct_eigenvalue`` — the latent consumer, now realized.**
  :func:`~orpheus.numerics.eigenvalue.direct_eigenvalue`'s dense
  ``np.linalg.solve(A, F)`` **is** this operator's action written as free
  functions — the latent consumer that made the promotion finishable
  (taxonomy §5: "no current consumer" often means "consumer not yet
  wired"). The engine itself stays **ndarray-pure** — its closed-form
  verification is the point — so it was never going to be the *production*
  spelling. That spelling landed at **task #138 (step 5b)**: the homogeneous
  solver now composes ``K = MatrixInverseOperator(loss) @ production`` and
  eigendecomposes ``K.as_matrix()``, making it the **first production
  consumer** of the class. It constructs the matrix inverse *explicitly*
  rather than via the structure-keyed ``loss.inverse()`` — which, reading the
  sum :math:`C - K_{\rm iso}` (invertible leading term), would return the
  iterative :class:`~orpheus.numerics.green_operator.GreenOperator` splitting
  — the direct-realization strategy encoded as a type. With the operator
  spelling live, ``direct_eigenvalue`` has **zero production consumers**,
  retained as the ``(A, F)``-posed sibling engine and the RQI test oracle.
  See :ref:`direct-eigensolve-solve`.
* **CP :math:`[P]` — the next production method in line.** The
  collision-probability matrix is dense **by construction**
  (:doc:`collision_probability`, §14b), so a CP ``.inverse()`` realized as a
  ``MatrixInverseOperator`` is the next production consumer after homogeneous.

This is the **same consumer ruling** ``GreenOperator`` shipped under at step
4: the type is correct and pinned, and its consumers arrive incrementally —
homogeneous first (step 5b), CP next. The *factory dispatch* that would route
``.inverse()`` to this type automatically is the remaining seam; explicit
construction by a size-aware consumer is the honest interim (and, for
homogeneous, the deliberate permanent choice).


Verification
------------

``MatrixInverseOperator`` is a **direct** inverse, so — unlike the iterative
:class:`~orpheus.numerics.green_operator.GreenOperator` — its correctness is
asserted at **machine·cond** grain against a **closed-form** dense reference
(hand-built ``Diagonal`` / ``Permutation`` matrices; the
structurally-independent
:meth:`~orpheus.transport.operators.isotropic_scattering.IsotropicScattering.dense_per_material`
*storage* transpose for the energy leaves), never at a driver tolerance. It
shares Green's honest framing in one respect: it is not an eigenvalue solver
and is not source-driven, so **no gate here makes an eigenvalue claim on an
MMS reference**; the homogeneous :math:`k_\infty` that rides through the
retired-loop relocation is anchored on its landed **closed-form** SymPy
value (``test_kinf_exact``, 1e-12), the analytical pillar.

.. list-table:: Step-5 verification gates (all ``@pytest.mark.foundation``)
   :header-rows: 1
   :widths: 42 58

   * - Gate file
     - What it pins
   * - ``tests/numerics/test_matrix_inverse_operator.py``
     - The base ``as_matrix`` L0 (exact vs hand-built + the storage-oracle
       cross-check; the **C-order column convention on a non-symmetric op**;
       rectangular-honesty; :math:`\equiv` the retired ``_as_dense`` loop),
       the ``ValueError`` / ``MatrixTooLarge`` boundary (class-discriminated,
       with an **at-threshold** designed-green control), and the
       ``MatrixInverseOperator`` invariants (M-materialise + M-direct at
       machine·cond with the **Green driver-tol contrast**, seed
       bit-identity, the back-half anchor, the non-square / singular /
       too-large guards with a positive constructs-cleanly control, and the
       :math:`(-S_{\rm ao})+D` **witness**).
   * - ``test_inverse_universal.py`` /
       ``test_operator_capability_predicates.py`` /
       ``test_operators_apply_typed.py``
     - The 4th sibling's participation in the universal family: the
       direct-construction registry row, object faithfulness, and the static
       ``assert_type`` / ``SupportsSeededApply`` conformance (a kwarg-less
       ``apply`` override is a CLI-pyright
       ``reportIncompatibleMethodOverride``).

The config discipline follows the family: the column convention is pinned on
a **non-symmetric** operator (a symmetric one is blind to a transpose,
Mode 6), and every value gate anchors on the **hand-built** reference, not on
``np.linalg.solve(A.as_matrix(), ·)`` alone (which would be self-referential).
The teeth are a **14-mutation bank** verified under ``-O`` (2026-07-02): each
mutation reddens a *named* gate — the ``as_matrix`` transpose / ravel /
size-gate / resolve family, the ``MatrixInverseOperator`` LU-transpose /
seed-consume / forward-``as_matrix`` / non-square-guard / kwarg-drop /
structural-guard-added family, the two ``OperatorProduct``-closure mutations,
and the homogeneous-relocation divergence — beside the designed-green
controls (the at-threshold materialization, the positive constructor, the
ignored seed). The pyright ratchet held exactly at 148; the homogeneous
:math:`k_\infty` / flux stayed byte-identical; the Sphinx build stayed
``-W`` clean.


.. _operator-algebra-assembly-axis:

The assembly axis — structural sparse emission (stencil-assembly 2b)
====================================================================

Step 5 established ``as_matrix`` as the functor **out** of the operator
category, :math:`\mathrm{Op}\to\mathrm{Mat}`, realised by apply-to-basis
probing (:eq:`matrix-functor-out`), and deferred a *sparse* override with
its consumer. Stencil-assembly 2b lands that override — as a full
**assembly axis** parallel to the inverse and adjoint axes. The
per-method cell mathematics (the SN closure walk, the diffusion FD
stencil) is developed in :doc:`loss_representations`
(:ref:`loss-rep-three-modes`) and :doc:`diffusion_1d`; this section is
the *operator-algebra* view: how emission threads through the composers,
how ``as_matrix`` delegates, and why the retained probing pathway is
kept as a permanent oracle.

Two realizations of the Mat-functor
-----------------------------------

The Mat-functor now has **two** realizations, and the split is the same
total-versus-partial distinction that separates ``as_matrix`` (a total
functor) from ``inverse()`` (partial — only where an inverse exists;
:ref:`capability-set-semantics`):

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - Realization
     - Domain
     - Cost
   * - **apply-to-basis probing** (``_as_matrix_by_probing``)
     - **total** — every finite-dimensional operator has a matrix
     - dense, :math:`n` operator applications
   * - **structural emission** (:meth:`~orpheus.numerics.operator.SupportsAssembly.assemble`)
     - **partial** — only where a leaf declares a stencil
     - sparse, :math:`O(\mathrm{nnz})` scatter

Probing is the honest fallback for any operator (it reads only
``apply``); emission is the efficient path for a structured operator
that knows its own :math:`(\text{row},\text{col},\text{value})`
footprint. The emitted matrix lands in
:class:`~orpheus.numerics.assembled_operator.SparseAssembledOperator` — a
thin :mod:`scipy.sparse` wrapper conforming to
:class:`~orpheus.numerics.operator.LinearOperator`, **not** a new
COO-builder type with its own algebra (that would twin the operator
algebra one layer down — every law restated on triplet buffers). scipy's
own ``COO → CSR`` conversion performs the FEM duplicate-summing, so an
emitter may scatter per-cell / per-face contributions freely and the
carrier assembles them. The functor is **closed** on that carrier: an
assembled operator's ``assemble()`` is itself.

The three-layer assembly surface
--------------------------------

Assembly is the third **per-axis three-layer surface**, minted exactly
like the inverse and adjoint axes (:ref:`capability-set-semantics`):

#. a **predicate** —
   :attr:`~orpheus.numerics.operator.LinearOperator.is_assemblable`
   (base default ``False``; the runtime, instance-accurate truth,
   recursive on composites) — the successor idea to the retired
   capability tags, reading structure not a registry;
#. a **narrowing target** —
   :class:`~orpheus.numerics.operator.SupportsAssembly` (a Protocol
   *extending* ``LinearOperator``, declaring ``assemble() ->
   SparseAssembledOperator``), deliberately **not** ``runtime_checkable``
   (an ``isinstance`` would match every ``OperatorSum``, which defines
   ``assemble`` class-uniformly even when a summand cannot emit);
#. a **checked bridge** —
   :func:`~orpheus.numerics.operator.assemblable` (a PEP-647
   ``TypeGuard``): the one construct that turns the runtime predicate
   into the static permission to call ``assemble()``.

The refusal is the assembly-axis sibling of
:class:`~orpheus.numerics.operator.NotInvertible` (inverse axis) and
:class:`~orpheus.numerics.operator.MissingAdjoint` (adjoint axis):
:class:`~orpheus.numerics.operator.MissingAssembly`, a ``TypeError``
raised **eagerly** by the composer ``assemble()`` bodies when an operand
cannot emit. An operator without a stencil simply does **not** declare
``assemble`` (misuse is a *static* error, method absence over an
advertising flag), and the probing ``as_matrix`` remains its total
Mat-functor. A space-anonymous leaf (a bare-ndarray multiplier / iso
operator with no composite layout) reports ``is_assemblable = False``
honestly — there is no global DOF numbering to emit into.

The composer homomorphism laws
------------------------------

Emission is an **additive-monoidal functor**, so a composite assembles
by recursion through its operands' emissions — no re-walk of the
stencils. Each composer carries one homomorphism law in its
``assemble()`` body, and the matching ``is_assemblable`` predicate is the
conjunction over its legs:

.. math::

   [A+B] = [A] + [B], \qquad
   [A\,B] = [A]\,[B], \qquad
   [\alpha L] = \alpha\,[L],

realised respectively by the carrier's own CSR **addition**
(:class:`~orpheus.numerics.operator.OperatorSum`), CSR **matmul**
(:class:`~orpheus.numerics.operator.OperatorProduct`), and scalar
**multiply** (:class:`~orpheus.numerics.operator.ScaledOperator`). A sum
or product is assemblable iff **both** legs are; a scaled operator iff
its operand is; and the eager :class:`~orpheus.numerics.operator.MissingAssembly`
guard-narrows the legs (Design C) before any operand call. The
:class:`~orpheus.numerics.operator.TensorProductOperator` law would be
:math:`[A\otimes B] = [A]\otimes[B]` (a :func:`scipy.sparse.kron`), but
it is **deferred with no consumer** — the diffusion loss and the SN
:math:`L(+C)` trees contain no tensor-product leaf, so building it now
would be a primitive with no product (Cardinal Rule 2).

R2 — ``as_matrix`` delegates to densified assembly
--------------------------------------------------

The ruling **R2** wires the two realizations together without changing
``as_matrix``'s contract. ``as_matrix`` keeps its dense semantics — same
basis-shape resolution, same :class:`~orpheus.numerics.operator.MatrixTooLarge`
size gate — and then, when the operator is
:func:`~orpheus.numerics.operator.assemblable`, **delegates** the
densification to ``self.assemble().as_matrix(...)`` (the assembled
matrix's own ``.toarray()``, with the column dimension checked against
the resolved basis shape) instead of running :math:`n` probing
applications; otherwise it falls back to the probing loop. Because **no
operator was assemblable until an emitter landed**, the delegation is a
no-op for every pre-2b call site — bit-safety by construction. The
:class:`~orpheus.numerics.flat_operator.FlattenedOperator` is transparent
to the axis: its ``is_assemblable`` / ``assemble()`` are its inner
operator's, since the typed operator's emission is already in the flat
layout.

The anti-tautology discipline
-----------------------------

Once ``as_matrix`` delegates to assembly, a naive
"``as_matrix`` :math:`\equiv` ``assemble().to_dense()``" gate compares
assembly **with itself** — vacuous. The defence is a retained,
**separately named** pathway: the probing loop lives on as
``_as_matrix_by_probing`` (not inlined in ``as_matrix``), precisely so a
gate can *force* the probing realization on an assemblable operator and
compare it against the structural emission. This is the
fuller-view-oracle discipline (the same reason the rolling-window sweep
keeps its full-field oracle): the retained relinquished view is the
verification pathway that pins the optimized path, and **an assembly bug
must never be able to hide inside its own densification**.

The G3 gate is exactly **one permanent pin per family** — the diffusion
loss (``test_g3_probed_equals_assembled_pin``) and one DD slab block
(``test_g3_dd_slab_probed_column_pin``) — forcing
``_as_matrix_by_probing`` against ``assemble().as_matrix()``. One pin per
family suffices: the equivalence is a structural property of the
delegation, not a per-instance numerical coincidence.

The first production consumer — the diffusion resolvent
-------------------------------------------------------

The assembly axis is not a primitive-without-a-product: its first
production consumer is the **diffusion resolvent**. The exact inner
solve is ``MatrixInverseOperator(FlattenedOperator(A, template))`` — one
eager LU at construction, one back-substitution per outer iteration
(:doc:`diffusion_1d`). Because every diffusion loss leaf now emits,
``A.as_matrix()`` routes through the R2 delegation, so
:class:`~orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator`
LU-factors the **assembled** matrix automatically — no consumer-side
change, and the whole existing diffusion suite (its k / trace / balance
gates) becomes the assembled path's regression net.

That "automatically" is the exact place a coverage gate can go vacuous
(vv-principles **Mode 11** — a green twin that never executes the
rewired line), so it is pinned by a **sentinel** rather than trusted:
``test_resolvent_materializes_through_assembly`` monkeypatches
``_as_matrix_by_probing`` to a counter and asserts that constructing the
resolvent fires **no probe at the composite flat dimension** — the
delegation to assembly genuinely executed. (Probes at the tiny law
dimension :math:`n_g` **are** expected and allowed: the boundary emitter
extracts each realized law's :math:`n_g\times n_g` block *through* the
law's own ``apply`` — that in-emitter probing is the one-source
discipline, not a delegation escape.) A green equivalence gate proves the
values; only the sentinel proves the consumer is on the new path. The
diffusion probed-versus-assembled densification measured a max
:math:`|\Delta| = 0.0` (bit-identical) on the heterogeneous fixture.


.. _coupled-block-operator:

The coupled block operator — the ψ½ ray as System B
===================================================

The curvilinear-S\ :sub:`N` within-group problem is the operator-surface
taxonomy's capstone: a **2×2 coupled block operator** over two systems.
It is the block-level generalisation of the three-layer surface
(:ref:`capability-set-semantics`) — ``apply`` becomes the block
matvec, ``assemble`` scatters each block at its DOF offset, and ``solve``
runs structure-keyed block substitution — and it is where the
curvilinear starting-direction flux :math:`\psi_{1/2}` finally earns a
first-class home. The campaign that built it (``coupled_block_operator``,
GH #280 the walk unification / #282 the direct ψ½ seed) replaced a
*fused* implementation — the ψ½ seed hand-rolled inside the sweep and
inside the model-generic scattering gain — with an explicit block system
in which every coupling is a named operator and a wrong pairing is
**unconstructable**. This section documents the posing, the four blocks,
the N-general machinery, the one production spelling, the block solve,
the convergence certificate, and the swap law. The starting-direction
physics (the pole as a straight characteristic, the M1/M2/M3 metric
distinction, R12a presence) lives in
:ref:`sn-282-direct-starting-direction-solve` in :doc:`discrete_ordinates`;
here we document the **algebra** the physics is posed in.

Why ψ½ is a system, not a kwarg — the two-point radial BVP
----------------------------------------------------------

At the two closed :math:`\mu = \pm 1` rays of a curvilinear μ-level the
angular-redistribution coefficient :math:`1-\mu^2` vanishes
(:math:`\alpha_{1/2}=0`, Hébert 3.423), so the streaming–collision
balance for :math:`\psi_{1/2}` **decouples** from the α-cascade and
reduces to a plain diamond-difference recurrence in radius —
:math:`\mu\,\partial_r + \sigma_t`, a *straight characteristic* (Hébert
§3.9.4). This ODE is not a scalar boundary datum the bulk sweep can
consume as a kwarg: it is a **two-point boundary-value problem** in its
own right, carrying **two** boundary conditions —

* **r = R Dirichlet** — the outer-face inflow corner (:math:`\mu = -1`)
  is *given* data: 0 for vacuum, the reflected outflow corner for
  reflective;
* **r = 0 pole continuation** — :math:`\psi_{1/2}^{+}(0) =
  \psi_{1/2}^{-}(0)`, the inward leg's exit face is the outward leg's
  entry face.

A quantity closed by its own two-point BVP is a **system**, not a
parameter. Posing it as one lets the algebra answer "what invariant
tests it?" (the march is the exact inverse of its own recurrence),
"how is it coupled?" (two named off-diagonal operators), and "how is it
solved?" (block substitution) — the four operator-algebra questions the
fused kwarg spelling could not even ask.

The within-group augmented system is therefore

.. math::
   :label: coupled-block-2x2

   \begin{bmatrix} A_{AA} & A_{AB} \\ A_{BA} & A_{BB} \end{bmatrix}
   \begin{bmatrix} \psi_A \\ \psi_B \end{bmatrix}
   \;=\;
   \begin{bmatrix} q_A \\ q_B \end{bmatrix},
   \qquad A_{ij} : V_j \to W_i ,

.. vv-status: coupled-block-2x2 documented

with **System A** the transport bulk ⊕ trace (the angular-flux
:class:`~orpheus.transport.full_field.FullField`, ``BulkField`` bulk
angular flux ⊕ its spatial ``BoundaryField`` trace) and **System B** the
ψ½ radial-characteristic ray (the ``RadialCharacteristicField`` cells at
each radial cell, carrying its two BC data). This is not a new tensor
type — it is the composite biproduct re-partitioned,
:math:`\mathrm{Mat}_2(\mathrm{Mat}_2(\mathcal C)) \cong
\mathrm{Mat}_4(\mathcal C)`. The two-system membership is carried by the
:class:`~orpheus.numerics.operator.SystemRole` axis (below), orthogonal to
the bulk↔boundary :class:`~orpheus.numerics.operator.BlockRole` that
refines System A.

The four blocks
---------------

Each block is a named production operator with its own code home,
invariants, and adjoint. The diagonal blocks are the systems'
self-operators; the off-diagonals are the couplings. All four live in
:mod:`orpheus.sn.operators.radial_characteristic` (System A's self-block
alone has no explicit object — it is the driver-level composite).

.. list-table:: The 2×2 blocks and their realisations
   :header-rows: 1
   :widths: 8 30 34 28

   * - Block
     - Object
     - What it is
     - Invariants
   * - :math:`A_{AA}`
     - the driver composite :math:`L + C - S - B_a`
     - System A's self-block: streaming + collision − bulk scattering
       gain − trace boundary gain
     - carries :class:`BlockRole` ``FULL`` by the join lattice;
       ``solve`` = the WDD sweep
   * - :math:`A_{AB}`
     - :class:`~orpheus.sn.operators.radial_characteristic.RadialCharacteristicSeeding`
     - ray → bulk seed injection (the Morel–Montry angular seed) —
       **cell-local angular**, no spatial cell-cell coupling
     - :math:`\sigma`-**independent** (a type fact); ``is_adjointable``,
       ``is_invertible`` ``False`` (rectangular)
   * - :math:`A_{BA}`
     - :class:`~orpheus.sn.operators.radial_characteristic.RadialCharacteristicEmission`
     - bulk → ray coupling, :math:`\mathrm{Fold}\circ K\circ\int d\mu`
       (the emission the within-group gain lags)
     - ``is_adjointable`` (its transpose IS the S-adjoint pullback);
       ``is_invertible`` ``False``
   * - :math:`A_{BB}`
     - :class:`~orpheus.sn.operators.radial_characteristic.RadialCharacteristicOperator`
     - the banded radial DD march :math:`\mu\,\partial_r + \sigma_t`
     - ``is_invertible`` = ``is_adjointable`` = ``True``; the direct
       Carlson march IS :math:`A_{BB}^{-1}`

**A naming gotcha to spell out.** ``A_BB`` (the class
:class:`~orpheus.sn.operators.radial_characteristic.RadialCharacteristicOperator`)
is the **bare** radial march :math:`\mu\,\partial_r + \sigma_t`. System
B's self-block in the *loss* grid is :math:`A_{BB} - B_b`, the march minus
the ray-corner boundary gain :math:`B_b`
(:class:`~orpheus.sn.operators.boundary.RadialCharacteristicBoundaryOperator`)
— **exactly parallel to System A**, whose self-block is
:math:`A_{AA} = L + C - S - B_a`. Both systems factor as
``(transport operator) − (bulk/reaction gains) − (boundary gain)``, and
both boundary gains :math:`B_a`, :math:`B_b` are lagged (they live in
:math:`N`, below), each reflecting through its **own** trace so the
:math:`-B` term arrives as given data (:ref:`bc-extraction`).

The self-block :math:`A_{BB}` — the march is its own inverse
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Because the recurrence is triangular in radius (:math:`\rho = 0` — the
#284 forward-substitution certificate, measured :math:`A_{ss} = 5.0`
self-coupling, :math:`A_{sb} = 0` exact), the two-leg Carlson march
(:func:`~orpheus.sn.spatial.psi_half_angle_seed.carlson_inward_sweep_from_source`,
Hébert 3.434–3.435) is the **exact** direct inverse :math:`A_{BB}^{-1}` —
no iteration, no previous-iterate seed. :meth:`RadialCharacteristicOperator.solve
<orpheus.sn.operators.radial_characteristic.RadialCharacteristicOperator.solve>`
marches the inward :math:`\mu=-1` leg from the r=R inflow corner to the
pole, then rides the *same* engine on reversed cell data (orientation
carried by the data, never a flag) from the pole-continued face out to
the r=R outflow corner. All four action surfaces are realised —
:meth:`apply` / :meth:`apply_transpose` (the forward
:math:`A_{BB}^{\mathsf T}`), :meth:`solve` / :meth:`solve_transpose`
(:math:`(A_{BB}^{-1})^{\mathsf T}`) — so ``is_invertible`` and
``is_adjointable`` both read ``True`` and the involution web
``inverse().solve == apply`` closes. The ``apply ∘ solve`` outflow-corner
defect closes to ``0.0`` bit-exactly; the cell round-trip is
principled-equivalent at ~FP ULP (the forward's :math:`2/\Delta r` and
the march's :math:`\Delta r\,\sigma + 2` reassociate). Since campaign
step 4e the production ``(L+C)`` walk routes its ray solve **through**
this operator — the former in-sweep inline march is retired, so the
two-leg orchestration lives in exactly one place (Cardinal Rule 2).

The seed block :math:`A_{AB}` — a σ-independent cell-local angular coupling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

At each radial cell :math:`i`, the ray value :math:`\psi_{1/2}(i)` seeds
the Morel–Montry angular recurrence (Hébert §3.9.4)

.. math::
   :label: coupled-ab-seed

   \psi_{m+1/2,\,i}
     = \frac{\psi_{m,\,i} - (1-\tau_m)\,\psi_{m-1/2,\,i}}{\tau_m},
   \qquad \psi_{-1/2,\,i} \equiv \psi_{1/2}(i),

run over ordinates :math:`m` at a **fixed** cell :math:`i`. The upstream
half-flux enters that cell's balance as the angular numerator
:math:`(\Delta A/w)\,c_{\rm in}\,\psi_{m-1/2,\,i}`. So the seed at cell
:math:`i` couples ONLY to the bulk ordinates at the *same* cell — there
is **no spatial cell-cell coupling** (that is :math:`A_{BB}`'s job). This
is what separates :math:`A_{AB}` from :math:`A_{BB}`: :math:`A_{AB}` is
cell-local angular and realises both directions HERE as thin wraps of the
single-sourced Morel–Montry closure (:meth:`apply` zeroes the bulk to
isolate its block by linearity; :meth:`apply_transpose` runs the local
gather :math:`\bar n = -\bar o/V` then reverses the recurrence to the seed
cotangent ``seed_cells_bar`` — the Euclidean transpose
:math:`A_{AB}^{\mathsf T}`). It needs **no** :math:`\sigma_t`: with the
bulk zeroed the collision/streaming terms drop out and only the
σ-independent angular numerator survives, so the coupling is a pure
function of the mesh geometry and quadrature. That σ-independence is a
**type fact** — the constructor takes only the mesh.

.. vv-status: coupled-ab-seed documented

The emission block :math:`A_{BA}` — the Schur fold ``Fold ∘ K ∘ integrate``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The bulk-within-group flux emits an isotropic source that seeds the ray.
:math:`A_{BA}`
(:class:`~orpheus.sn.operators.radial_characteristic.RadialCharacteristicEmission`)
composes three factors,

.. math::
   :label: coupled-ba-emission

   A_{BA} \;=\; \underbrace{\mathrm{Fold}}_{\text{Reconstruction}}
            \;\circ\; \underbrace{K}_{\text{emission kernel}}
            \;\circ\; \underbrace{\textstyle\int\! d\mu}_{\text{integrate}} ,

the angular integral of the bulk flux to :math:`\phi_0`, the operator's
isotropic emission kernel :math:`K`, and the reconstruction of that
:math:`\ell = 0` moment at the closed rays. The kernel is a **dependency
injection**: :math:`A_{BA}` is generic over any ``ndarray → ndarray``
emitter carrying ``apply``/``apply_transpose``. The production
instantiation passes the scattering
:attr:`~orpheus.transport.operators.scattering.ScatteringOperator.isotropic_kernel`
(:math:`K_{\rm iso} = \Sigma_{s0} + 2\Sigma_{2n}`) — the **same shared
object** the bulk scattering gain uses, so the emission is single-sourced
(one shared kernel call, never a twin re-implementation of
:math:`\Sigma_{s0}^{\mathsf T}\phi`).

.. vv-status: coupled-ba-emission documented

The **fold factor** :math:`\mathrm{Fold}`
(:class:`~orpheus.sn.operators.radial_characteristic.RadialCharacteristicReconstruction`)
is the 1-D angular reconstruction (Hébert Eq. 3.432) *sampled* at the
closed rays :math:`\mu = \pm 1`:

.. math::
   :label: coupled-ba-fold

   \bar q_{1/2}(\mu = \pm 1)
     \;=\; \sum_\ell \frac{2\ell + 1}{2}\,q_\ell\,(\pm 1)^\ell ,

the same :math:`\tfrac{2\ell+1}{2}\,(\pm 1)^\ell` weights the angular
frame reconstructs with (:math:`P_\ell(\pm 1) = (\pm 1)^\ell`), evaluated
at the rays rather than the quadrature nodes. Both the forward fold and
its Euclidean transpose route through the ONE math kernel
:func:`~orpheus.numerics.spaces.radial_characteristic_space.fold_moments_to_radial_characteristic`
(Cardinal Rule 2 — the :math:`P_1(-1) = -1` sign is spelled once). The
transpose is the exact **broadcast/sum adjoint pair**: the forward folds
one moment source onto every carried level and both signs (a broadcast);
the transpose *sums* the per-level, per-sign ray cotangents back into
moment space. Corners stay zero — the fold is volumetric; the
inflow-corner datum is :math:`B_b`'s job.

.. vv-status: coupled-ba-fold documented

**The LIFT — S/F pure-bulk, the driver composes the emission.** Before
campaign step 4c the model-generic scattering/fission composites
hand-rolled the ψ½ seed inside their own ``apply`` — a curvilinear-SN
augmentation welded into a model-generic gain. The **lift** (the Wave-O
#208 pattern that separated :math:`B` from :math:`S`) reversed that: it
made :math:`S`/:math:`F` **pure bulk** and posed the coupling as a
first-class operator the driver lags as a separate gain. The
within-group scattering gain rides :math:`(S,\ A_{BA},\ B_a)`, and
:math:`A_{BA}`'s transpose is exactly the ``w · K_isoᵀ(Reconstructionᵀ
χ_seed)`` bulk pullback the S-adjoint used to carry inline — moved HERE,
single-sourced, so ``S.apply_transpose`` is pure bulk and
:math:`(L+C-S-A_{BA}-B).H` reconstructs the monolithic adjoint. Because
the forward keeps :math:`K_{\rm iso}` while the adjoint rides the fold
transpose, the reciprocity gate is a genuine cross-check of two
structurally-different representations of one operator, not a tautology.

.. note::

   **Fission does NOT flow through** :math:`A_{BA}` **(HAZARD 5).** The
   emission is kernel-generic and accepts the fission dyad
   :math:`\chi \otimes \nu\Sigma_f` as a smoke-verified second kernel, but
   fission's *production* ray seed rides the outer :math:`q_{\rm ext}`
   seam as a **direct** :math:`\mathrm{Fold}`: fission is the eigenvalue
   outer source, so its :math:`K\circ\int d\mu` is already pre-computed as
   the fission source, and routing it through the full
   :math:`\mathrm{Fold}\circ K\circ\int d\mu` would **double-apply**
   :math:`K\circ\int`. Within-group fission is zero (it enters as
   :math:`q_{\rm ext}` per the eigenvalue-outer / within-group split), so
   :math:`A_{BA}`'s production shape is the *scatter* coupling; the
   genericity keeps the emission a clean dependency injection, not a claim
   that fission wires through it.

The N-general block machinery
-----------------------------

The 2×2 ψ½ system is **instance #1** of a semantics-agnostic N-system
machinery in :mod:`orpheus.numerics.coupled_system` — nothing there knows
transport, rays, or meshes. A coupled system is N sub-systems solved
together: the state is the direct sum of the systems' fields
(:class:`~orpheus.numerics.coupled_system.CoupledField`), the space is
their product (:class:`~orpheus.numerics.coupled_system.CoupledSpace`),
and the operator is the N×N grid of blocks
(:class:`~orpheus.numerics.coupled_system.CoupledOperator`) with the block
matvec

.. math::
   :label: coupled-block-matvec

   y_i \;=\; \sum_j A_{ij}\, x_j

as the *only* spelling. A missing block (``None``) **is** the zero map,
structurally — block existence doubles as coupling sparsity, and no
zero-padding arithmetic ever runs.

.. vv-status: coupled-block-matvec documented

**Why a typed grid over a padded** ``OperatorSum`` **(the rejected
route).** The alternative — a flat sum of same-space operators, each
padding the blocks it does not touch with present-zeros — was **rejected**
(campaign re-scope, 2026-07-10): padding keeps *wrong multiplications
representable*. A padded block accepts any composite; nothing at the type
level stops a ray operator from receiving a bulk field or an emission
from landing in the trace slot, and ``system_role`` tags are runtime
metadata, not type constraints. The honest object is a typed block
operator over a typed block vector, where every term of
:eq:`coupled-block-matvec` type-checks per block and a wrong pairing is
**unconstructable** (``coding-elegance`` Pattern 1 ∘ Pattern 4 — the
algebra is the syntax).

**The three block-level modes** mirror the three-layer operator surface,
one level up: ``apply`` (the block matvec, Krylov action), ``assemble``
(each block's sparse emission scattered at its ``(row_i, col_j)`` DOF
offset into one flat matrix via :func:`scipy.sparse.block_array`, closing
over :class:`~orpheus.numerics.assembled_operator.SparseAssembledOperator`
— the same axis as :ref:`operator-algebra-assembly-axis`), and ``solve``
(the structure-keyed direct solve, below). The offsets ARE the block
structure — :attr:`CoupledSpace.system_slices
<orpheus.numerics.coupled_system.CoupledSpace.system_slices>` is the
scoped local→global DOF map ``assemble`` needs, the same layout
:meth:`CoupledField.to_flat` packs. Because :class:`CoupledField` itself
satisfies the ravellable ``to_flat``/``from_flat`` protocol, every
``restart = n_dof = template.to_flat().size`` sizing site tracks the
coupled dimension automatically (the ERR-053 GMRES-truncation family is
closed by conformance, not per-site edits).

**The Hilbert adjoint comes free — Mode-12 closure by construction.**
:class:`CoupledOperator` implements only the Euclidean
:meth:`apply_transpose` — the transposed grid :math:`(A^{\mathsf T})_{ji}
= (A_{ij})^{\mathsf T}` — and carries a :class:`CoupledSpace` whose metric
methods dispatch member-wise. The metric adjoint :math:`A.H = G^{+}
A^{\mathsf T} G` is then realised ONCE by the existing
:class:`~orpheus.numerics.operator._AdjointOperator` wrapper. No adjoint
code lives in the block machinery — which is exactly what keeps the block
adjoint Mode-12-closed: a hand-rolled "Euclidean block ``.H``" that skips
the metric conjugation is the ERR-067 reopening (``vv-principles``
Mode 12). The block ``.H`` inheriting the member metrics is the campaign's
highest-value verification row.

The two-system role lattice — :class:`SystemRole`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:class:`~orpheus.numerics.operator.SystemRole` tags which of the two
systems an operator's action maps between —
:attr:`A` (within System A), :attr:`B` (within System B: the self-block
:math:`A_{BB}` and the ray boundary :math:`B_b`), or :attr:`COUPLED` (an
off-diagonal, or the assembled 2×2: :math:`A_{AB}`, :math:`A_{BA}`, and
the grid). Reading each role as the *set* of systems its action touches
(:math:`A = \{A\}`, :math:`B = \{B\}`, :math:`\text{COUPLED} = \{A,B\}`),
a sum touches the union: "same role stays, anything different becomes
COUPLED". This is the two-system analogue of the bulk↔boundary
:class:`~orpheus.numerics.operator.BlockRole` join, with ``COUPLED``
playing the top-of-lattice role ``FULL`` plays there — a deliberate twin
kept while only two role axes exist (the collapse trigger is a third
parallel axis, e.g. a DSA/multiphysics role). Model-generic leaves
(``C``/``S``/``F``, every diffusion/CP/MoC operator) leave ``system_role``
at ``None`` — the conservative "not part of the ψ½ augmentation" reading —
so a System-A composite must be **stamped** ``SystemRole.A`` explicitly at
its composition site (the model-generic members' honest ``None`` would
otherwise poison the join to ``None``).

The one production spelling — ``build_within_group_system``
-----------------------------------------------------------

The joint system has exactly ONE construction site:
:func:`~orpheus.sn.coupled_system.build_within_group_system`, which
returns the frozen :class:`~orpheus.sn.coupled_system.WithinGroupSystem`
record — the loss grid together with its **named regular splitting**
:math:`A = M - N` (Hackbusch 2016 §11), all four members built from the
SAME piece objects (one ``L+C``, one ``S``, one ``B_a``, one ``B_b``, …).
This builder is what retired the former ``_within_group_triple`` /
``_lagged_gains`` construction pair. The grid and its
:class:`CoupledSpace` are emitted **together**, aligned by construction
(RULING P1: a mismatched operator/space pairing is unconstructable).

The **loss grid** carries the loss-sign convention, and the two
off-diagonals differ in sign — a trap worth spelling out:

.. math::
   :label: coupled-loss-grid

   A \;=\;
   \begin{bmatrix}
     L + C - S - B_a & +\,\text{Seeding} \\[2pt]
     -\,\text{Emission} & A_{BB} - B_b
   \end{bmatrix}

The :math:`(A,B)` seed is **positive** (its ``apply`` already emits the
seed's term of the fused bulk residual row — the loss sign is *internal*
to the operator, matching the in-sweep placement it wraps); the
:math:`(B,A)` emission is **negated** (a
:class:`~orpheus.numerics.operator.ScaledOperator` — the emission is a
gain, so the ray equation :math:`(A_{BB} - B_b)\psi_B -
\text{Emission}\,\psi_A = q_B` carries the minus).

.. vv-status: coupled-loss-grid documented

**The splitting.** The SI/Krylov drivers do not consume the loss grid
raw — they consume its regular splitting :math:`A = M - N`,

.. math::
   :label: coupled-mn-splitting

   M = \begin{bmatrix} L+C & +\,\text{Seeding} \\ \mathbf 0 & A_{BB}
       \end{bmatrix},
   \qquad
   N = \begin{bmatrix} S + B_a & \mathbf 0 \\ +\,\text{Emission} & B_b
       \end{bmatrix},

with :math:`M` the **sweepable part** inverted every step and :math:`N`
the lagged coupling gains (all signs positive — gains on the rhs
:math:`\text{rhs} = q + N\psi`; the loss grid's :math:`-\text{Emission}`,
:math:`-B_b` minus signs are the :math:`M - N` complement, not a
contradiction). Note :math:`M`'s :math:`(B,B)` entry is the **bare** march
:math:`A_{BB}` while :math:`N`'s is :math:`B_b`, so
:math:`A(B,B) = A_{BB} - B_b` recovers the loss self-block. The drivers
iterate :math:`\psi \leftarrow M^{-1}(q + N\psi)` (SI) or GMRES on
:math:`(M-N)\psi = q`.

.. vv-status: coupled-mn-splitting documented

**Presence is structural (R12a).** System B exists **only** on a
seed-carrying mesh (the 1-D sphere): the mesh's ray spaces are ``None``
:math:`\iff` non-carrying, and the ``RadialCharacteristic*`` block
constructors *refuse* seedless meshes. So a seed-carrying mesh builds the
2×2; a seedless mesh (slab / cylinder / Cartesian) builds the 1×1
:math:`[[A_{AA}]]` and the splitting degrades to the bare :math:`(L+C,\
(S, B_a))` the seedless driver paths consume zero-touch. "Applying a
System-B block on a non-carrying mesh" is not a runtime branch — it is an
object that **does not exist**. Since the B.2d eviction ``FullField`` is a
pure 2-block composite (:math:`\psi_A` cannot carry ray state at all), the
old dead-slot double-count hazard dissolved structurally: the coupled flat
dimension is the honest two-system sum, no dead padding.

The block-triangular solve and the materialised-LU EXTRACT
----------------------------------------------------------

On a carrying mesh :math:`M` is an **honest upper-triangular**
:class:`CoupledOperator` grid, so its ``solve`` is matrix-free block
back-substitution — and the substitution order is exactly the curvilinear
sweep order:

.. math::
   :label: coupled-block-substitution

   \psi_B \;=\; A_{BB}^{-1}\,q_B ,
   \qquad
   \psi_A \;=\; (L+C)^{-1}\bigl(q_A - \text{Seeding}\,\psi_B\bigr).

System B's march runs first (its exact direct inverse — one pass, no
seed), then System A's bulk sweep runs on its **ray-decoupled** channel
:math:`q_A - \text{Seeding}\,\psi_B`. This is one body for all four
orientation × transpose combinations: the transpose of a triangular grid
is triangular the other way (:math:`(A^{\mathsf T})_{ij} =
(A_{ji})^{\mathsf T}`), so the same substitution runs with the visit
order flipped, guarded per block by
:class:`~orpheus.numerics.operator.MissingAdjoint`. The inverse **object**
is :class:`~orpheus.numerics.coupled_system.CoupledSubstitutionOperator`
(the taxonomy §12 wrap-delegate — its ``apply`` delegates to the grid's
``solve``, ``initial_guess`` accepted and dropped since a direct
substitution has nothing to seed).

.. vv-status: coupled-block-substitution documented

**The materialise/LU EXTRACT (step 5a).** The *full* 2×2 loss grid is not
triangular (:math:`A_{BA} \neq 0`), so it takes the second route:
:meth:`CoupledOperator.inverse` materialises the assembled matrix and
LU-factors it via
:class:`~orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator`
(one eager LU, then back-substitution per solve). This EXTRACT is
**principled-equivalent, not bit-identical** to the production
splitting-iteration — the matrix reduction tree differs — and a naive
extraction returns O(1) garbage (the sweep treats inflow/seed rows as
*given data*, so the row-contract must be preserved). It is the oracle the
swap-law gates ride; production solves stay the splitting iteration on the
record's ``resolvent``/``gains``. The **iterative** splitting solve
(block-Jacobi / block-Gauss-Seidel over :math:`A = M - N`) deliberately
stays with the drivers — convergence is spectral
(:math:`\rho(M^{-1}N) < 1`), never a structural capability.

Why the :math:`(B,A)` block is ``None`` in :math:`M` — the Schur/lag argument
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ray↔bulk coupling splits cleanly by *which direction iterates*. The
seed :math:`A_{AB}` sits **inside** :math:`M` (the :math:`(A,B)` slot),
because within :math:`(L+C)` the ray→bulk coupling is block-triangular and
solved **directly** in one pass — measured :math:`A_{sb} = 0` exact (no
bulk→seed feedback inside :math:`L+C`), :math:`A_{bs} = 7.5` (the M-M
seed feed). The emission :math:`A_{BA}` sits **outside** the resolvent, in
:math:`N` — it is the bulk *scattering* gain's ray arm, and ONLY the
scattering coupling iterates: measured :math:`S_{sb} = 0.183` (bulk→ray
source), :math:`S_{bs} = 0` exact (ψ½ carries zero moment weight), with
outer spectral radius :math:`\rho(M^{-1}N) = 0.371 < c` (Adams–Larsen).
Placing :math:`A_{BA}` in :math:`M` would fold a lagged gain into the
one-pass resolvent — it belongs on the rhs, lagged. So the :math:`(B,A)`
slot of :math:`M` is structurally zero and the resolvent stays
block-**triangular** (direct); the coupling that genuinely iterates rides
:math:`N`. "Deciding where the block goes IS operator algebra": the same
four posed operators would let :math:`A_{BA}` sit folded, explicit, or in
a DSA preconditioner — a composition choice the machinery supports.

The ρ-honest stop and the lag-death certificate
-----------------------------------------------

Splitting the emission into a lagged gain creates a verification hazard:
an iteration whose fixed point does **not** solve the equation (a stale or
lagged block inside :math:`M` — the #282 class) can still report
"converged" to any :math:`\lVert\Delta\psi\rVert`-family stop. Step 5c
closes this in two pieces.

**The running stop is the ρ-honest equation residual** (a deliberate
re-interpretation of ``tol`` — :class:`~orpheus.numerics.iteration.SourceIteration`):

.. math::
   :label: coupled-free-identity-residual

   r_{n} \;=\; A\,\psi_{n} - q_{\rm ext}
         \;=\; \mathrm{rhs}_{n-1} - \mathrm{rhs}_{n}
         \;=\; \textstyle\sum_i g_i\,(\psi_{n-1} - \psi_{n}),
   \qquad
   {\rm res}_n = \frac{\lVert r_n\rVert_2}
                      {\max(\lVert q_{\rm ext}\rVert_2,\ 10^{-30})} .

The **free-identity** spelling :math:`\mathrm{rhs}_{n-1} -
\mathrm{rhs}_{n}` (retained from the loop's own bookkeeping — zero
marginal cost, checked before the next apply) is exact when the step
operator is an exact inverse of :math:`M`: then :math:`M\psi_n =
\mathrm{rhs}_{n-1}` and :math:`A\psi_n - q = N(\psi_{n-1} - \psi_n)`. The
residual is preferred over :math:`\lVert\Delta\psi\rVert` because the
increment understates the true error by :math:`1/(1-\rho)` (Banach) AND is
blind to the lag-death class; the residual claim is
contraction-rate-independent.

.. vv-status: coupled-free-identity-residual documented

**The certificate closes the exact-**:math:`M` **assumption.** The
free identity itself *assumes* exact :math:`M` — the hole the driver-level
**convergence certificate**
:class:`~orpheus.sn.solver.ConvergenceCertificateError` closes. Once, at a
*claimed* exit, :func:`~orpheus.sn.solver.evaluate_residual` measures the
true :math:`r = A\psi - q` through a real forward apply — the only
measurement an in-:math:`M` lag cannot fool — and raises if the defect
exceeds :math:`10\times\text{tol}`. It is wired on every full-angular arm
(the coupled sphere, the seedless un-windowed SI, both Krylov paths). The
headline gate is a stale-zero-:math:`\psi_B` in-:math:`M` lag mutation: the
running stop reports convergent while the certificate raises
``match="lag-death"`` — with control legs on both sides, the classifier's
asymmetry proof (#282 measured the cold-residual defect at ~5e5). The
typed residual carries System B as its **own** member
(:class:`~orpheus.transport.residuals.radial_characteristic_interior_residual.RadialCharacteristicInteriorResidual`
⊕ its boundary sibling), so a System-A-only residual cannot silently drop
a wrong seed row — the Mode-12 (b) closure. This is an additive diagnostic
+ certificate, NOT in the convergence path (see :ref:`affine-typed-residual`).

The swap law on the grid
------------------------

The inverse-as-operator taxonomy's swap law :math:`A.H.\text{inverse}()
\equiv A.\text{inverse}().H` extends to the grid: the joint adjoint IS the
grid's transposed substitution. :class:`CoupledSubstitutionOperator`
advertises ``is_adjointable`` iff every coupling block transposes AND every
diagonal carries the direct ``solve_transpose`` verb (the #280 two-factor
discipline per block), and its :meth:`apply_transpose` delegates to
:meth:`CoupledOperator.solve_transpose` (the transposed block
substitution). The step-6 collapse sharpened the single-system predicate
to match: :attr:`SweepOperator.is_adjointable
<orpheus.sn.operators.sweep_operator.SweepOperator.is_adjointable>` is now
the two-factor ``isinstance(inner, InvertibleOperator) and
inner.is_adjointable`` — the B.2d carrying-mesh *third* factor retired,
because with no ψ½ legs anywhere a bare :math:`(L+C)` on a carrying mesh IS
unambiguously the ray-decoupled :math:`(A,A)` block (the type says what you
hold). The joint adjoint's one home is the grid's transposed substitution;
there is no adjoint SI driver in production, so the :math:`A_{BA}`/:math:`A_{AB}`
pullbacks are reached only by the ``.H`` reciprocity gates — which is why
they are realised now and pinned by **nonzero-seed-cotangent** gates (a
present-zero seed would hide a lost pullback).

The one-channel collapse (step 6) — the walk IS the :math:`(A,A)` block
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

With the four blocks posed and the solve routed through the grid, the
walk's fused ψ½ joint-leg channel became dead production code and was
**retired** (step 6, net −803 lines). The two forward/transpose kwarg
pairs, the three-function presence-guard family, the walk's in-solve
System-B engines, and the operator-free ``transport_sweep`` wrapper are
all gone. Every walk surface is now the ray-decoupled :math:`(A,A)`
diagonal block, on **every** mesh, with no seed parameters: the matvec
substitutes a zero seed into the Morel–Montry thread (bit-identical to the
retired dead-slot arithmetic — the closure reads zeros), and the transpose
**discards** the thread cotangent (a fixed zero input's cotangent
propagates nowhere; the :math:`\text{Seeding}^{\mathsf T}` pullback is the
explicit grid block). The eigenvalue finalize re-routes through the SAME
:func:`~orpheus.sn.coupled_system.build_within_group_system` ``.resolvent``
every driver consumes. The mesh remains the single authority on presence;
what changed is that nothing *checks* against it anymore — the type system
carries the biconditional. The narrative of the walk's own view of this
collapse lives in :doc:`loss_representations`.

Numerical evidence
------------------

The block architecture is a **principled-equivalent** re-association of
the fused implementation — same operator algebra, different reduction
tree — so the sphere re-baselines are FP-grain, sphere-only, and the
seedless geometries stay bit-identical (the strongest leak tripwire).

.. list-table:: Measured equivalence (the fused → block-sum re-association class)
   :header-rows: 1
   :widths: 40 30 30

   * - Check
     - Result
     - What it pins
   * - :math:`A_{BB}` forward-substitution
       :math:`\lVert\text{solve}\circ\text{apply}(\psi)-\psi\rVert`
     - :math:`3.5\text{e-}16`
     - the march is the exact direct inverse (#284, :math:`\rho=0`)
   * - block-triangular certificate :math:`A_{sb}`
     - :math:`0` exact
     - no bulk→seed feedback inside :math:`L+C` (direct one-pass)
   * - EXTRACT vs dense-LU of assembled :math:`(L+C)`
     - :math:`5.5\text{e-}16`
     - the materialise/LU route is the principled oracle
   * - scatter LIFT vs the old monolithic ``S.apply``
     - :math:`1.2\text{e-}15`
     - :math:`S_{\rm bulk} + A_{BA}` ≡ the pre-lift composite
   * - finalize equivalence — slab / cylinder, every member
     - **bitwise**
     - seedless geometries are untouched (leak tripwire)
   * - finalize equivalence — sphere :math:`k`, scalar, both ψ½ members
     - **bitwise**
     - :math:`A_{BB}`'s march ≡ the retired in-solve march exactly
   * - finalize equivalence — sphere bulk + trace movers
     - rel-max :math:`8.4\text{e-}16`
     - the M-substitution vs fused-interleave re-association only

The outer spectral radius :math:`\rho(M^{-1}N) = 0.371` and the measured
block algebra (:math:`A_{bs} = 7.5`, :math:`A_{ss} = 5.0`, :math:`S_{sb} =
0.183`, :math:`S_{bs} = 0`) confirm the direct/iterative split the API
reflects. The full within-group wall (tests/sn + tests/numerics, not-slow
serial) is 3080/0 through the step-6 collapse.

The DSA seam — issue #2 as a future consumer
--------------------------------------------

The coupled-block machinery is not a primitive-without-a-product beyond
the ψ½ instance: **consistent DSA (Issue #2)** is its next consumer. A
diffusion-synthetic-acceleration step poses the SN transport system and a
low-order diffusion correction as a **coupled system** — an N-system
:class:`CoupledOperator` whose diagonal blocks are the transport loss and
the diffusion loss operator :math:`A_{\rm diff} = L + C - S - B` (which
already exists, :mod:`orpheus.diffusion.operators`, #290 P4) and whose
off-diagonals are the restriction :math:`R` (transport residual → coarse)
and prolongation :math:`P` (coarse correction → fine) operators. Those
:math:`R`/:math:`P` operators **do not exist yet** — they are the future
``SystemRole`` members Issue #2 will add; the seam is *named*, not built.
When they land, DSA consumes the same three modes (the block matvec, the
``assemble`` route the diffusion resolvent already rides via
:ref:`operator-algebra-assembly-axis`, and the block ``solve``), and it
consumes the transport residual :func:`~orpheus.sn.solver.evaluate_residual`
already types — DSA computes that residual, the diffusion solve corrects
it, and the correction :math:`\to 0` at convergence (so the accelerator is
correctness-safe by construction). The Wave-T DSA row
(:ref:`wave-t-tensor-network`) tracks the same seam from the
loss-leaf side.


.. _trace-spaces-doc:

Trace spaces — :math:`\Gamma_-` and :math:`\Gamma_+`
====================================================

Boundary conditions act on the **directional half** of the
transport equation's boundary trace. Per Grand Report v3 §5.3 and
§16A.5, the trace splits into two pieces by the sign of
:math:`\Omega \cdot \hat n`:

.. math::
   :label: trace-half-decomposition

   \Gamma_- \;=\; \{(\mathbf{r}, \Omega) \in \partial\Omega \times S^d
                  : \Omega \cdot \hat n(\mathbf{r}) < 0\},
   \qquad
   \Gamma_+ \;=\; \{(\mathbf{r}, \Omega) : \Omega \cdot \hat n > 0\}.

.. vv-status: trace-half-decomposition documented

The two halves are represented by typed
:class:`~orpheus.numerics.space.FunctionSpace` subclasses
:class:`~orpheus.numerics.spaces.angular_trace_space.InflowTraceSpace` and
:class:`~orpheus.numerics.spaces.angular_trace_space.OutflowTraceSpace`, which
carry a **per-face directional mask**:

.. math::
   :label: per-face-inflow-mask

   \mathrm{inflow\_mask}[f, n]
   \;=\;
   \bigl(\Omega_n \cdot \hat n_f < -\epsilon\bigr).

.. vv-status: per-face-inflow-mask documented

The mask has shape ``(n_faces, n_ordinates)`` boolean; the
tangential band :math:`|\Omega_n \cdot \hat n_f| \leq \epsilon`
(default ``ε = 1e-12``) is in neither half.

Three consumers read the mask today:

* The SN realizer's vacuum branch consumes
  ``inflow_indices_for_face(face)`` to construct an
  :class:`~orpheus.numerics.operator.IncomingOrdinateMaskTensor`
  with the per-face inflow indices.
* The universal invariant
  :meth:`~orpheus.geometry.boundary.BoundaryTraceLaw.assert_source_lives_on_incoming_trace`
  uses the inflow mask to validate a
  :class:`~orpheus.geometry.boundary.InflowSourceSpec` (ERR-047).
* The SN curvilinear sweep (1-D spherical / cylindrical) consumes
  the same realizer-routed mask as the slab and 2-D Cartesian
  paths (Issue #188 + #176, closed 2026-05-11).

Construction goes through the classmethod factory
:meth:`AngularTraceSpace.from_quadrature_and_layout
<orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.from_quadrature_and_layout>`,
which (since Issue #225 / C5.3) is **geometry-blind**: it builds the
per-face :math:`\Omega\cdot\hat n` masks from the angular quadrature's
direction-cosine arrays and a :class:`~orpheus.numerics.face_layout.FaceLayout`
whose ``"{axis}{min|max}"`` face names imply axis-aligned outward
normals — no spatial mesh is consulted. Every :class:`Mesh1D` coord
system (``CARTESIAN`` / ``SPHERICAL`` / ``CYLINDRICAL``) shares the
same ``("xmin", "xmax")`` radial-axis face structure —
:class:`GaussLegendre1D` is the shared quadrature, with ``mu_x`` as the
direction cosine along that axis — and 2-D / 3-D Cartesian add the
``y`` / ``z`` faces from the same convention. The 2-D cylindrical
(axisymmetric :math:`(r, z)`) case never reaches the factory: such a
:class:`Mesh2D` cannot become an :class:`SNMesh` (no 2-D cylindrical SN
sweep exists), so the refusal lives at the :class:`SNMesh` construction
surface, not the trace factory. See :ref:`sn-c5-geometry-blind-trace`.

The two trace spaces and the
:class:`~orpheus.geometry.boundary.InflowSourceSpec` Protocol
together close the §16A.1 affine boundary form
:math:`\gamma_- \psi = R\,G\,\gamma_+ \psi + q` documented in
detail at :ref:`affine-bc-form`.


.. _eigenvalue-posing:

Eigenvalue posing and the power-iteration algorithm
===================================================

The operator leaves :math:`L, C, S, F, B` documented above are the raw
material; this section is the **assembly instruction** that turns them
into a criticality eigenproblem and runs it. The architecture is
**layered into four tiers**, and the layering is the load-bearing
design decision: it is what makes the :math:`\alpha`-eigenvalue,
adjoint, and transient problems land later as *pure additions* (new
posing data) rather than new solver engines. The corrected layering
was confirmed by an independent structural analysis (the
``cross-domain-attacker`` ``eigenvalue_posing_layering_frames`` and
``power_iteration_vs_keigenvalue_morphism`` memos, 2026-06-04) and
realized in commits ``650032e`` / ``7603c8e`` (2026-06-05).

.. admonition:: Key Facts (eigenvalue posing)
   :class: tip

   - **Standard form:** the generalized eigenproblem
     :math:`A_{\rm loss}\,\psi = \lambda\,M\,\psi`. Its power-method
     realization is the dominant eigenpair of the **resolvent**
     :math:`A_{\rm loss}^{-1} M`.
   - **Krein–Rutman / Perron–Frobenius:** for a compact, positive
     :math:`A_{\rm loss}^{-1} M` the fundamental mode is the *unique*
     non-negative eigenvector and the dominant eigenvalue is real and
     positive — the only physically meaningful steady state. All
     higher harmonics change sign in space.
   - **k-eigenvalue (LIVE):** :math:`A_{\rm loss} = L+C-S-B`,
     :math:`M = F`, :math:`k = \mu`. The dominant eigenvalue of
     :math:`A_{\rm loss}^{-1} F` is :math:`k_{\rm eff}`.
   - **Four layers:** operator leaves (method-specific) → problem
     posing (bifurcated 2a/2b) → resolvent :math:`A_{\rm loss}^{-1}`
     (method-specific) → solution algorithm (general over the standard
     form).
   - **The invariant:** the Layer-4 algorithm sees ONLY a
     normalized-source fixed-point procedure — *apply* :math:`M`,
     *solve* :math:`A_{\rm loss}^{-1}`, *normalize*, *estimate* the
     dominant :math:`\mu`. It never touches the method's operators or
     sweeps.
   - :func:`~orpheus.numerics.eigenvalue.power_iteration` is the
     **canonical** Layer-4 algorithm (the *more general* layer);
     :class:`~orpheus.numerics.iteration.KEigenvalue` is one Layer-2b
     implementer that delegates its loop to it. **One loop.**


The standard form and its resolvent
------------------------------------

Discretizing the steady-state transport (or diffusion) equation
produces a balance between **loss** (streaming + collision − in-group
scattering − boundary in-scatter) and **production** (fission). Group
every loss term into a single operator :math:`A_{\rm loss}` and every
production term into a single **eigen-operator** :math:`M`. The
criticality condition — that a self-sustaining flux distribution
exists when production is scaled by :math:`1/\lambda` — is the
**generalized eigenproblem**

.. (vv-status rationale) The generalized-eigenproblem standard form.
   The verifiable claim — that the dominant eigenvalue of
   A_loss⁻¹ M equals k_eff for the k-posing — is anchored against the
   homogeneous closed-form algebra-of-record
   (orpheus.derivations.common.eigenvalue.kinf_and_spectrum_homogeneous,
   k = λ_max(A⁻¹F)); the production iteration is the L1/L2 chain of the
   five family solvers.
.. vv-status: eigen-standard-form documented

.. math::
   :label: eigen-standard-form

   A_{\rm loss}\,\psi \;=\; \lambda\,M\,\psi .

Inverting the loss operator turns the generalized eigenproblem into a
standard one for the **resolvent operator**
:math:`K_{\rm pm} \equiv A_{\rm loss}^{-1} M`:

.. math::
   :label: eigen-resolvent

   K_{\rm pm}\,\psi \;=\; \tfrac{1}{\lambda}\,\psi ,
   \qquad
   K_{\rm pm} \;\equiv\; A_{\rm loss}^{-1} M .

This is the form every power-method realization actually iterates: one
outer step is *apply* :math:`M`, then *invert* :math:`A_{\rm loss}`
(the fixed-source solve), then *renormalize*, then *estimate*
:math:`\lambda`. The reason the dominant eigenpair is the one the
iteration converges to — and the reason it is the *physical* one — is
the **Krein–Rutman theorem** (the infinite-dimensional Perron–Frobenius
statement): for a compact, positive :math:`K_{\rm pm}`,

* the dominant eigenvalue :math:`\rho(K_{\rm pm})` is real, positive,
  and simple;
* its eigenvector is the *unique* eigenvector with no sign changes —
  i.e. a physically realizable non-negative flux distribution;
* every higher harmonic changes sign in space and is therefore not a
  steady reactor state.

Power iteration converges to exactly this fundamental mode, at a rate
governed by the **dominance ratio** :math:`|\lambda_1/\lambda_0| =
|k_1/k_0|` (Trefethen & Bau 1997, §27). The dominant K and
:math:`\alpha` eigenvalues are *extreme* eigenvalues of
:math:`A_{\rm loss}^{-1} M`, reachable by plain power iteration;
**shift-invert** :math:`(A_{\rm loss} - \sigma M)^{-1}` is the strict
generalization needed only for *interior* eigenvalues (higher
harmonics, FEAST/Arnoldi–Schur), and is a documented future seam, not
a present need.

.. note::

   The standard form :eq:`eigen-standard-form` is the discrete twin of
   the algebra-of-record in
   :func:`~orpheus.derivations.common.eigenvalue.kinf_and_spectrum_homogeneous`,
   which solves :math:`k = \lambda_{\max}(\mathbf{A}^{-1}\mathbf{F})`
   with :math:`\mathbf{A} = \mathrm{diag}(\Sigma_t) - (\Sigma_s +
   2\Sigma_2)^{T}` and :math:`\mathbf{F} = \chi \otimes \nu\Sigma_f`
   for the homogeneous infinite medium. That closed-form reference is
   the **structurally-independent** ground (a closed-form analytical
   pillar, exact in the homogeneous limit; reducing to the 1-group
   :math:`k = \nu\Sigma_f/\Sigma_a`) against which the production
   power iteration's converged eigenvalue is verified. The agreement
   of two ORPHEUS solvers is cross-implementation agreement, not
   correctness evidence; the closed-form ratio is the anchor.


The four layers
---------------

.. list-table:: The four-tier eigenvalue architecture
   :header-rows: 1
   :widths: 6 22 16 56

   * - Layer
     - Name
     - Specificity
     - What lives here
   * - 1
     - Operator leaves
     - method-specific
     - :math:`L, C, S, F, B` (+ future :math:`T = 1/v`). The
       block-diagonal :math:`G`-metric (codomain inner product of
       :math:`L`'s composite :math:`V_{\rm bulk}\oplus V_{\rm trace}`
       space: bulk :math:`V\,w_n` :math:`\oplus` trace
       :math:`|\Omega\cdot\hat n|\,w_n`) lives HERE and is reused by
       every posing (step O.2b R5, :ref:`bc-extraction-g-adjoint`).
   * - 2
     - Problem posing
     - **bifurcated**
     - **2a** (method-agnostic): role assignment + the :math:`\mu \to`
       physical-eigenvalue map — which leaves play :math:`A_{\rm loss}`
       vs :math:`M`, and how :math:`\mu` maps to :math:`k` / :math:`\alpha`.
       **2b** (method-specific): how the method assembles and inverts
       the concrete :math:`A_{\rm loss}` object.
   * - 3
     - Resolvent :math:`A_{\rm loss}^{-1}`
     - method-specific
     - the fixed-source inner solve. SN:
       :class:`~orpheus.numerics.iteration.SourceIteration` /
       :class:`~orpheus.numerics.iteration.KrylovAcceleration`. CP:
       BiCGSTAB. Diffusion: direct FD inverse (exact LU of the fused
       :math:`A`, no inner iteration — #290). Inverts *whatever*
       :math:`A_{\rm loss}` the posing produced; independent of problem
       type.
   * - 4
     - Solution algorithm
     - general over the standard form
     - eigenvalue-finders (power iteration | full-spectrum Arnoldi /
       Krylov–Schur | shift-invert / FEAST) over
       :math:`(A_{\rm loss}^{-1}, M)`; time-integrators (transient)
       over :math:`(A_{\rm loss}, T, q(t))`.

**Why posing bifurcates (2a vs 2b).** The first-draft architecture
treated posing as wholly method-agnostic — "just arrange the leaves."
That is false: the *role assignment* is agnostic, but the *loss-operator
realization* is method-specific. SN realises its loss operator as the
invertible resolvent :math:`L + C` (the WDD sweep) plus the lagged
coupling gains :math:`S` (bulk scattering) and :math:`B` (boundary
reflection), handed to the **variadic** within-group driver as
:math:`(L+C,\,S,\,B)` (:ref:`bc-extraction-variadic-driver`). CP has
**no** :math:`(A, S, F)` split at all — its
:meth:`solve_fixed_source <orpheus.numerics.eigenvalue.EigenvalueSolver.solve_fixed_source>`
is one BiCGSTAB on a *monolithic* collision-probability matrix; the
factor :math:`(A-S)^{-1}` does not exist as a separable object.
Splitting the posing into

* **2a — role assignment + μ-map** (pure data: a posing-table row), and
* **2b — loss-operator realization** (the method's concrete assembly),

lets :math:`2a \circ 2b \circ 3 \circ 4` compose cleanly across every
family. The key consequence:
:class:`~orpheus.numerics.iteration.KEigenvalue(A, S, F)` is the
operator-triple **2b realization** — NOT a problem-type layer. Treating
the operator triple as a "problem type" was the conflation the
bifurcation removes.

**The variadic driver IS the posing/resolvent boundary made explicit.**
The Layer-3 SN resolvent
(:class:`~orpheus.numerics.iteration.SourceIteration` /
:class:`~orpheus.numerics.iteration.KrylovAcceleration`) is now
**problem-type-agnostic**: it consumes
:math:`\text{Driver}(A_{\rm resolvent},\,*\text{gains})` and never asks
which gain plays which posing role. *Which* leaves are gains is the
2a decision — for the SN k-row the gains are exactly the
:math:`A_{\rm loss}` couplings :math:`S` and :math:`B` (fission
:math:`F` is the eigen-operator :math:`M`, not a within-group gain; it
enters as :math:`q_{\rm ext}`). The retired fixed :math:`(A, S, F)`
triple had baked a 2a role distinction into the Layer-3 resolvent,
where it does not belong — the variadic generalisation pushes the
distinction back up to the posing layer (:ref:`bc-extraction-variadic-driver`).

**The invariant (Layer-4 sees only a fixed point).** Layer 4 consumes
the method-agnostic
:class:`~orpheus.numerics.eigenvalue.EigenvalueSolver` boundary — a
normalized-source fixed-point procedure with exactly these moves:
build the eigen-source from the current flux, solve the resolvent,
renormalize to unit production rate, estimate the dominant eigenvalue,
test convergence. It never sees SN sweeps, CP matrices,
:class:`~orpheus.transport.timed_full_field.TimedFullField` carriers,
or angular state — those live at Layers 1–3, *below* the boundary. This
is the abstraction that makes a new problem type a posing-row addition
rather than an engine rewrite.


The posing table
-----------------

Each row assigns the leaves to roles, gives the eigen-operator, and
gives the :math:`\mu \to` physical-eigenvalue map. **The k-row is LIVE;
the rest are documented future seams** — recorded here so a future
session lands them as pure additions.

.. list-table:: Eigenvalue posing rows (2a — method-agnostic)
   :header-rows: 1
   :widths: 16 30 14 26 14

   * - Problem type
     - :math:`A_{\rm loss}`
     - eigen-operator :math:`M`
     - :math:`\mu \to` physical map
     - Status
   * - **k-eigenvalue**
     - :math:`L + C - S - B`
     - :math:`F`
     - :math:`k = \mu`
     - **LIVE**
   * - :math:`\alpha`-eigenvalue
     - :math:`L + C - S - F - B`
     - :math:`T = 1/v`
     - :math:`\alpha = -1/\mu`
     - future seam
   * - adjoint
     - :math:`A_{\rm loss}^{\dagger}` (daggered row)
     - :math:`M^{\dagger}`
     - same :math:`\lambda` as the forward row
     - future seam
   * - transient
     - :math:`A_{\rm loss} + T/\Delta t` (per implicit step)
     - — (source-driven, not eigen)
     - — (time integration)
     - future seam

**The k-row, in full.** For k-eigenvalue criticality the fission
production is the eigen-operator and *everything else* is loss. The
resolvent's dominant eigenvalue is :math:`k_{\rm eff}` directly
(:math:`k = \mu`):

.. math::
   :label: eigen-k-posing

   (L + C - S - B)\,\psi \;=\; \tfrac{1}{k}\,F\,\psi
   \qquad\Longleftrightarrow\qquad
   \bigl[(L+C-S-B)^{-1} F\bigr]\,\psi \;=\; k\,\psi .

This is exactly :eq:`operator-eigenvalue` with the boundary in-scatter
:math:`B` made explicit (Wave O step O.4a.2 promoted :math:`B` to a
first-class sibling leaf; see :ref:`bc-extraction`). In production the
within-group loss :math:`L+C-S-B` is realised honestly: :math:`S` and
:math:`B` are two separate coupling gains handed to the variadic driver
(Wave O step O.2a — :ref:`bc-extraction-variadic-driver`), so the
matvec is :math:`(L+C).\text{apply} - S.\text{apply} - B.\text{apply}`.
The transitional :math:`S + B` driver fold is retired.

**The α-row (future seam).** The :math:`\alpha`-eigenvalue (the
time-eigenvalue, governing the asymptotic exponential time behaviour)
follows from the ansatz :math:`\psi(\mathbf r, \Omega, t) \propto
e^{\alpha t}\,\psi(\mathbf r, \Omega)`. Substituting into the
time-dependent transport equation, the time derivative
:math:`\frac1v \partial_t \psi` becomes :math:`\frac{\alpha}{v}\psi`,
so the steady balance reads

.. math::
   :label: eigen-alpha-derivation

   \tfrac{\alpha}{v}\,\psi + (L + C - S - F)\,\psi \;=\; 0
   \qquad\Longleftrightarrow\qquad
   (L + C - S - F - B)\,\psi \;=\; -\alpha\,T\,\psi ,
   \quad T \equiv \tfrac1v .

Matching to the standard form :eq:`eigen-standard-form` gives
:math:`A_{\rm loss} = L+C-S-F-B` (fission now joins the *loss* side,
because it is no longer the eigen-operator), :math:`M = T = 1/v`, and
:math:`\mu = -1/\alpha`, i.e. :math:`\alpha = -1/\mu`. The only new
machinery the :math:`\alpha`-row needs is a sixth leaf — a
:class:`~orpheus.numerics.operator.DiagonalOperator` realizing
:math:`T = 1/v` — joining :math:`L, C, S, F, B`. The posing, resolvent,
and algorithm layers are unchanged; this is the cleanest possible fit
and is why the layering was designed this way. **Not built** (only K
exists; *unify-after-two*).

**The adjoint row (future seam).** The adjoint eigenproblem
:math:`A_{\rm loss}^{\dagger}\,\psi^{\dagger} = \lambda\,M^{\dagger}\,
\psi^{\dagger}` is **just another posing row** whose role-operators are
the daggers of the forward leaves. The dagger is *free* from the
dagger-biproduct category already documented on this page (the ``.H``
adjoint propagates through ``+`` / ``&`` / ``@`` — see
:eq:`tensor-product-adjoint-distributivity`) and, crucially, from the
**composite metric-correct G-adjoint** that step O.2b R5 made concrete:
``op.H`` is the metric-correct :math:`A^{\dagger} = G^{-1} A^{\mathsf T}
G` over the :math:`V_{\rm bulk}\oplus V_{\rm trace}` carrier
(:ref:`bc-extraction-g-adjoint`), dense-oracle-verified to round-off.
Because forward and adjoint share the spectrum, :math:`\lambda` — and
therefore the :math:`\mu \to` physical map — is **unchanged**; only the
operators are daggered. The adjoint slots in at 2a with **zero new
engine machinery**: it is a row, not a layer, and its loss-operator
dagger is *already built and verified* by R5. (The first-draft
architecture's instinct to make adjoint a separate "mode" is the same
conflation the 2a/2b split removes.)

**The transient row (future seam).** Backward-Euler time stepping
:math:`(T/\Delta t + A_{\rm loss})\,\psi^{n+1} = (T/\Delta t)\,\psi^{n}
+ q^{n+1}` is a fixed-source solve with a **shifted** loss operator
:math:`A_{\rm loss} + T/\Delta t`. It **shares the Layer-3 resolvent**:
:meth:`solve_fixed_source <orpheus.numerics.eigenvalue.EigenvalueSolver.solve_fixed_source>`
inverts whatever loss operator the posing hands it, and the shifted
operator is still a streaming-plus-collision-like invertible object the
same sweep / BiCGSTAB handles. Transient therefore needs only (a) a
transient posing row and (b) a *time-integrator* Layer-4 sibling of
:func:`~orpheus.numerics.eigenvalue.power_iteration` that loops the
fixed-source solve in time and advances delayed-neutron precursors.
**No new resolvent, no new leaves** beyond the :math:`T` leaf the
:math:`\alpha`-row already introduces.


Why ``power_iteration`` is canonical, not deprecated
----------------------------------------------------

An earlier framing carried a :class:`DeprecationWarning` on
:func:`~orpheus.numerics.eigenvalue.power_iteration`, intending to
migrate all five solver families onto
:class:`~orpheus.numerics.iteration.KEigenvalue`. **The deprecation
arrow pointed the wrong way** and was removed in ``650032e`` /
``7603c8e``. The two are the **same fixed-point combinator** — one
power-method loop (the five-step *build source → solve resolvent →
renormalize → estimate → converge?* body) — instantiated at two
different layers:

* :func:`~orpheus.numerics.eigenvalue.power_iteration` exposes the
  inner resolvent **late**, behind the opaque
  :meth:`EigenvalueSolver.solve_fixed_source <orpheus.numerics.eigenvalue.EigenvalueSolver.solve_fixed_source>`
  Protocol method — a morphism the solver owns.
* :class:`~orpheus.numerics.iteration.KEigenvalue` binds the resolvent
  **early**, building it as :math:`(A-S)^{-1}` from the operator triple
  via an inner :class:`~orpheus.numerics.iteration.SourceIteration`.

The late-bound layer is **strictly more general**: it admits *both* the
sweep-posed operator-triple resolvent (SN, MoC — where the inner
:math:`(A-S)^{-1}` is a source iteration over the sweep :math:`A^{-1}`)
*and* the **monolithic-matrix resolvent** (CP, diffusion, homogeneous —
a single direct inverse). The early-bound layer can only express
methods whose resolvent factors as :math:`(A-S)^{-1}` from a triple via
that inner sweep — strictly narrower. Diffusion is the instructive case:
it *does* now carry an in-algebra :math:`(L, C, S, F)` family (#290;
:doc:`/theory/diffusion_1d`), yet it still belongs to the monolithic
camp, because it has **no sweep** — its resolvent is the explicit
:class:`~orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator`
of the *fused* :math:`A = L + C - S - B` (the scattering already
subtracted in), not an :math:`(A-S)^{-1}` iterated over a within-group
solve. Forcing CP or diffusion into the narrow layer would mean
manufacturing a fictitious sweep they do not have. Therefore the
**Protocol layer is canonical** and the **triple layer is a
specialization that adapts into it**.

.. list-table:: ``power_iteration`` vs ``KEigenvalue`` — same morphism, two layers
   :header-rows: 1
   :widths: 22 39 39

   * -
     - :func:`~orpheus.numerics.eigenvalue.power_iteration`
     - :class:`~orpheus.numerics.iteration.KEigenvalue`
   * - Layer
     - 4 (algorithm) over the 2-boundary
     - 2b (operator-triple posing realization)
   * - Resolvent binding
     - **late** — opaque ``solve_fixed_source``
     - **early** — :math:`(A-S)^{-1}` from the triple
   * - Admits
     - SN, MoC, CP, diffusion, homogeneous (any Protocol implementer)
     - SN / MoC only (needs an :math:`(A,S,F)` triple)
   * - The loop
     - **owns** the single power-iteration loop body
     - **delegates** to ``power_iteration``
   * - Role
     - the canonical engine
     - one implementer of the boundary

After the fix, the loop body lives in **one place**
(:func:`~orpheus.numerics.eigenvalue.power_iteration`).
:class:`~orpheus.numerics.iteration.KEigenvalue` realizes the
:class:`~orpheus.numerics.eigenvalue.EigenvalueSolver` boundary from its
triple — :meth:`compute_fission_source <orpheus.numerics.iteration.KEigenvalue.compute_fission_source>`
:math:`= F\psi/k`,
:meth:`solve_fixed_source <orpheus.numerics.iteration.KEigenvalue.solve_fixed_source>`
:math:`= (A-S)^{-1} q` via the warm-started inner
:class:`~orpheus.numerics.iteration.SourceIteration`, the **hardwired**
:math:`k`- and production-estimators
(:meth:`~orpheus.numerics.iteration.KEigenvalue.compute_keff` /
:meth:`~orpheus.numerics.iteration.KEigenvalue.compute_production_rate`;
the pre-R8 injection kwargs retired at #259 P1), and the
:math:`\ge 3`-iteration :math:`dk`/:math:`d\phi` convergence test — then
:meth:`solve <orpheus.numerics.iteration.KEigenvalue.solve>` simply
calls ``power_iteration(self, max_iter=self.max_outer)``. SN production
(:func:`~orpheus.sn.solver.solve_sn`), CP, diffusion, MoC, and
homogeneous all drive the same loop directly via the Protocol; the
``KEigenvalue`` adapter is for callers who *have* a natural
:math:`(A,S,F)` triple and want to skip writing a full solver class.

.. note::

   The "``KEigenvalue`` regresses :math:`P_\ell` (anisotropic
   scattering)" objection from the migration era dissolves under this
   framing: :class:`~orpheus.numerics.iteration.SourceIteration` is
   type-agnostic and **angular-capable** — it routes the RHS through
   the ravellable protocol, so a typed
   :class:`~orpheus.transport.operators.scattering.ScatteringOperator` acting on a
   :class:`~orpheus.transport.full_field.FullField` (the history-blind
   operator carrier; the driver feeds its
   :class:`~orpheus.transport.timed_full_field.TimedFullField` iterate, which
   reaches the arm via MRO and carries the full angular flux on its bulk)
   carries :math:`P_\ell` correctly. The observed regression was a property of
   an L1 *test adapter* that collapsed angular flux to scalar between
   outer iterations (dropping the angular moments), not of
   ``KEigenvalue``. The decisive — and sufficient — reason
   ``KEigenvalue`` cannot be the universal engine is the
   CP/diffusion/homogeneous **monolithic-resolvent (no-sweep)** fact
   alone.


The metric lives at the leaf, not the posing
---------------------------------------------

The :math:`G`-metric is the *codomain inner product* of :math:`L`'s
composite :math:`V_{\rm bulk}\oplus V_{\rm trace}` space (step O.2b R5,
:ref:`bc-extraction-g-adjoint`): a block-diagonal Hilbert metric with a
bulk phase-space block :math:`V_{\rm cell}\,w_n` and a partial-current
trace block :math:`|\Omega\cdot\hat n|\,w_n`. It is **intrinsic to the
streaming leaf** — :math:`L` carries its ``domain`` / ``codomain``
:class:`~orpheus.numerics.spaces.full_field_space.FullFieldSpace` with
the per-block ``inner_product_weights``, and the ``.H`` adjoint reads
them via the unchanged
:class:`~orpheus.numerics.operator._AdjointOperator` wrapper. It is
**NOT a posing-layer concern**: posing *arranges* leaves; the leaves
already *know* their metric through their composite space. The
:math:`G`-weighting is applied **once at the op level** — the
:class:`~orpheus.numerics.operator._AdjointOperator` wrapper reads
:math:`G` off the composite ``domain`` / ``codomain`` of the *summed*
operator and conjugates the whole sum
:math:`G^{-1}(\cdot)^{\mathsf T} G`, never re-applying it per leaf
(:ref:`bc-extraction-g-adjoint`). Since P4.5 W-D the bulk leaves
:math:`C, S, F` advertise the same composite ``full_field_space`` as
:math:`L` / :math:`B` (so the *forward* within-group guard validates),
but their own ``apply_transpose``, where defined, is the metric-blind
Euclidean transpose — the metric is layered on once by the sum's
adjoint wrapper, so a leaf carrying the composite domain is no
double-application risk. Consequently the adjoint posing row gets the
correct :math:`G`-weighted transpose for free — the dagger functor
applied to a composite that already carries the metric — which is
precisely why the adjoint row adds no new machinery.


Honest scope — what is agnostic today, and what is not
------------------------------------------------------

The architecture is realized **minimally**: only the k-row and
:func:`~orpheus.numerics.eigenvalue.power_iteration` exist. Two honest
caveats record where the present code stops short of the ideal so a
future session does not mistake intent for fact.

* **The Layer-4 loop is not yet *literally* K/α-agnostic.** Today the
  eigenvalue scaling lives in the K-specific
  :meth:`compute_keff <orpheus.numerics.eigenvalue.EigenvalueSolver.compute_keff>`
  (production/absorption ratio) and the :math:`/k` placement in
  :meth:`compute_fission_source <orpheus.numerics.eigenvalue.EigenvalueSolver.compute_fission_source>`.
  :func:`~orpheus.numerics.eigenvalue.power_iteration` is agnostic only
  by *delegating* the estimate to the problem's ``compute_keff``.
  Making the loop literally agnostic (relocating the scaling into the
  algorithm as a Rayleigh-quotient update on
  :math:`A_{\rm loss}^{-1} M`, adding an ``apply_loss`` method, and
  renaming the K-flavoured Protocol methods to ``eigen_operator`` /
  ``mu_to_eigenvalue``) touches all five families' Protocol surface and
  is **the first step of the α-wave**, snapshot-bit-identity-gated. It
  is deferred because only K exists (premature to unify;
  *unify-after-two*).

* **The full-spectrum / shift-invert seam is reserved, not built.**
  :attr:`KEigenvalue.eigenvalue_method <orpheus.numerics.iteration.KEigenvalue.eigenvalue_method>`
  selects the Layer-4 algorithm. Only ``"power"`` is implemented; any
  other value raises :class:`NotImplementedError` at construction time.
  Full-spectrum Arnoldi / Krylov–Schur and shift-invert / FEAST (for
  interior eigenvalues — higher spatial harmonics) slot in at this
  exact dispatch point, consuming the same
  :math:`(A_{\rm loss}^{-1}, M)` boundary.

.. warning::

   The :math:`\alpha`-eigenvalue, adjoint, transient, and full-spectrum
   rows are **documented future seams, not implemented features**.
   There is zero :math:`\alpha` / transient / Arnoldi / shift-invert
   scaffolding in production transport today (the
   :mod:`orpheus.kinetics` solver is 0-D point kinetics, not a
   deterministic-transport :math:`\alpha`/transient solver). The
   layering exists so each lands as a pure addition — a new posing row
   (α / adjoint), a new leaf (:math:`T = 1/v`), and at most a new
   Layer-4 sibling (the transient time-integrator) — never a rewrite of
   the engine, the resolvent, or the existing leaves.


Verification status
--------------------

The discriminating gate for the canonical-loop refactor is
``tests/numerics/test_iteration.py::test_keigenvalue_matches_solve_sn_2g_slab``:
it stays green after
:class:`~orpheus.numerics.iteration.KEigenvalue` delegates to
:func:`~orpheus.numerics.eigenvalue.power_iteration`. This is the
**same-morphism evidence** — if the two had been different algorithms,
routing ``KEigenvalue``'s loop through ``power_iteration`` would change
the converged answer; bit-stable agreement on a **2-group** slab (≥2
groups is mandatory — a 1-group eigenvalue is the flux-shape-independent
ratio :math:`k = \nu\Sigma_f/\Sigma_a` and detects no operator error)
confirms they are one combinator at two layers. Because the
:class:`~orpheus.numerics.eigenvalue.EigenvalueSolver` Protocol and the
five family solvers are untouched by the refactor, **every** family's
eigenvalue snapshot is trivially bit-identical across the change.

The converged *eigenvalue* (a solver-level claim) is verified against
the closed-form algebra-of-record
:func:`~orpheus.derivations.common.eigenvalue.kinf_and_spectrum_homogeneous`
(:math:`k = \lambda_{\max}(\mathbf{A}^{-1}\mathbf{F})`) for the
homogeneous reflective limit — a structurally-independent closed-form
pillar, not a code-to-code comparison.


Development history
===================

This is a reverse-chronological (latest first) changelog of the major
**architectural** milestones in the operator algebra. Iteration-rate
work, gate counts, and intermediate replans are deliberately omitted —
see the GitHub issues and the per-phase plan files for that granularity.
Entries marked *(in development)* live on an unmerged feature branch and
have no landed merge-to-``main`` hash yet; trust ``git`` over this table
for merge status.

.. list-table::
   :header-rows: 1
   :widths: 10 50 12 28

   * - When
     - Architectural milestone
     - Issue
     - Where
   * - in dev
       (2026-07-04)
     - **The assembly axis — structural sparse emission as the third
       Mat-functor realization** (stencil-assembly campaign, Phase 2b).
       ``as_matrix`` gains a second realization beside apply-to-basis
       probing: leaves that know their stencil emit a
       :class:`~orpheus.numerics.assembled_operator.SparseAssembledOperator`
       (a :mod:`scipy.sparse` serialization of the operator, not a new
       COO-builder algebra), and the axis carries the full three-layer
       surface minted like inverse/adjoint —
       :attr:`~orpheus.numerics.operator.LinearOperator.is_assemblable`
       predicate, :class:`~orpheus.numerics.operator.SupportsAssembly`
       narrowing Protocol, :func:`~orpheus.numerics.operator.assemblable`
       ``TypeGuard`` bridge, and the eager
       :class:`~orpheus.numerics.operator.MissingAssembly` refusal. The
       composers recurse through the homomorphism laws (Sum → ``+``,
       Product → ``@``, Scaled → scalar ``·``; TensorProduct → ``kron``
       deferred, no consumer), and ``as_matrix`` **delegates** to the
       densified emission when :func:`~orpheus.numerics.operator.assemblable`
       (ruling R2), with the probing loop retained as
       ``_as_matrix_by_probing`` — the fallback AND the anti-tautology
       oracle the probed≡assembled gates force. First production consumer:
       the diffusion resolvent
       ``MatrixInverseOperator(FlattenedOperator(A, template))`` now
       LU-factors the **assembled** matrix automatically (probed↔assembled
       measured bit-identical, max :math:`|\Delta| = 0`; a Mode-11 sentinel
       proves the delegation executes). See
       :ref:`operator-algebra-assembly-axis`.
     - #272
     - *(in development)*
       ``refactor/spatial-promotion-assembly``
   * - in dev
       (2026-06-25)
     - **The carrier grid recognised as a double category, and the
       operator type made two-parameter** (Frame-projection campaign,
       P4.5). The transport carriers are identified as the cells of a
       :math:`(\text{Representation} \times \text{Role})` **double
       category** — horizontal 1-morphisms are the representation-changing
       frame faces :math:`M`/:math:`R` (role-generic), vertical
       1-morphisms are the role-changing cross sections :math:`C`,
       :math:`\Lambda`, :math:`F` (representation-generic), and scattering
       :math:`S = \tfrac{1}{W}(R\circ\Lambda\circ M) =
       \texttt{frame.conjugate}(\Lambda)` is the **2-cell**, with the
       existing 0-ULP windowed-vs-full crosscheck recognised as its
       interchange-law coherence witness (:ref:`carrier-grid-double-category`).
       The operator Protocol is widened to the honest two-parameter
       :class:`~orpheus.numerics.operator.LinearOperator` ``Protocol[Domain,
       Codomain]`` (``apply(x: Domain) -> Codomain``; the PEP-696 default
       ``Codomain = Domain`` keeps ``[V] ≡ [V, V]`` for the endomorphic
       majority; requires-python raised to ``>=3.13``) — **a grid cell IS an
       operator's** ``(Domain, Codomain)`` (:eq:`carrier-grid-operator-typing`).
       The accompanying structural finding — that a fully-typed
       ``Carrier[Representation, Role]`` is impossible (Role-arithmetic and
       Representation-shape each force a class; a parameterized carrier
       breaks the runtime units gate via erasure), so the flat MI leaves
       are the unique normal form — is documented at
       :ref:`carrier-grid-flat-leaf-normal-form`. **Realisation status:**
       the two-parameter operator type and the double-category framing have
       **landed on the branch**, and W-C/W-F have **confined the composite
       operator carrier to the timeless** :class:`~orpheus.transport.full_field.FullField`:
       the ``apply`` dispatch arms register on ``FullField`` (a driver's
       :class:`~orpheus.transport.timed_full_field.TimedFullField` iterate
       reaches them via MRO), and W-F realigned the
       heteromorphic-``apply`` ``@overload`` stubs
       (:ref:`heteromorphic-apply-typing`) for scattering/fission from
       ``TimedFullField`` to ``FullField`` to match that registration
       (typing-only; runtime byte-identical). The ``@overload`` confessions
       themselves are **kept** — they are the honest per-carrier surface, not
       a wart to retire; their deeper *dissolution* into a thin ``match``
       router over single-sourced primitives is parked on #261 (the
       Pattern M-vs-``match`` spelling question, :ref:`heteromorphic-apply-typing`).
       The secondary-carrier ``ScalarFlux`` arm is likewise **kept** as the
       typed entry-point for the #205 cross-method scalar consumers; the bare-
       ``ndarray`` arm is K-eigenvalue-live. The deeper *secondary-carrier-arm*
       collapse couples to the C/F/S core relocation and CP / MoC carrier
       unification (#261).
     - #65 / #268 / #261
     - *(in development)*
       ``refactor/operator-inverse-algebra``
   * - 2026-06-08
     - **Flux states typed as an affine space; the iterate increment is a
       typed displacement.** ``flux − flux`` mints a
       :class:`~orpheus.transport.displacements._displacement.Displacement`,
       ``flux ⊕ displacement`` is the torsor update, and ``flux + flux`` is
       a :class:`TypeError` — the #201 dimensional gate becomes a *type*
       consequence (:ref:`affine-typed-field-algebra`). The Role axis of
       the carrier grid.
     - #208 / #201
     - ``main`` (Wave O step O.2)
   * - 2026-06-03
     - **The boundary law** :math:`B` **becomes a first-class sibling
       operator** and every operator's ``.apply`` output is typed as a
       *source/sink* (the bulk → :class:`~orpheus.transport.source_sinks.angular_source_sink.AngularSourceSink`,
       the boundary → :class:`~orpheus.transport.source_sinks.angular_boundary_source_sink.AngularBoundarySourceSink`),
       completing the operator-output role typing
       (:ref:`bc-extraction-operator-output-typing`).
     - #208
     - ``main`` (Wave O steps O.4a.2 / B.5.2)

