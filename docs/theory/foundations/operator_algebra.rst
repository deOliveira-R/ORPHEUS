.. _operator-algebra:

=============================
Operator Algebra Architecture
=============================

ORPHEUS' transport, eigenvalue, and Krylov solvers all act on a flux
distribution :math:`\psi` by composing a small set of linear operators,
each with a **distinct intrinsic mathematical type**:

- **collision** :math:`C = M[\sigma_t]` — a *multiplication operator*
  (diagonal; pointwise multiplication by the total cross section);
- **scattering** :math:`S = R\circ\Lambda\circ M` and **fission**
  :math:`F = |\chi\rangle\langle\nu\Sigma_f|` — *integral kernels*
  (nonlocal: scattering redistributes in angle, fission is the rank-1
  emission dyad);
- **streaming** :math:`L = \hat\Omega\cdot\nabla` and the **boundary
  law** :math:`B` — the leakage and its trace closure.

They compose into the **within-group transport operator**

.. (vv-status rationale) The governing within-group composition
   A = L+C−S−B — the loss operator. Definitional identity; the
   assembled composite is exercised by build_within_group_system and
   the fixed-source / eigenvalue suites, matching the sentineled
   operator-fixed-source / operator-eigenvalue siblings.
.. vv-status: operator-within-group-composition documented

.. math::
   :label: operator-within-group-composition

   A \;=\; L + C - S - B ,

posed either as an eigenvalue problem :math:`A\psi = \tfrac{1}{k}F\psi`
or a fixed-source problem :math:`A\psi = q`. Here :math:`A` is the
**loss operator** — removal and leakage net of the within-group gains —
against the fission **gain** :math:`F`. The boundary law :math:`B` is a
**first-class sibling** (not folded into :math:`L`), and :math:`F` sits
on the **right-hand side** — never inside :math:`A`. The invertible
sub-composite :math:`L+C` — streaming leakage plus total collision — has
the transport **sweep** as its exact inverse; :math:`S` and :math:`B`
are the within-group gains the outer iteration lags.

The :mod:`orpheus.numerics.operator` module installs these as a uniform
*matrix-free* algebra, so the eigenvalue, fixed-source, and
preconditioned-Krylov code consumes any method (S\ :sub:`N` / MoC / CP /
diffusion) without knowing which transport discretisation lives
underneath. This page is that algebra's reference: the intrinsic type of
each operator, the composition laws, the invertible :math:`L+C` and its
:term:`sweep`, and how a method extends the shared operators (S\ :sub:`N`
expands :math:`S` for anisotropy).

.. contents::
   :local:
   :depth: 2


Key Facts
=========

- **The realized boundary law** :math:`B` **is a first-class sibling
  operator** (Issue #208): the within-group transport
  operator is :math:`A = L + C - S - B` on the **two-block** transport
  state :math:`V = V_{\rm bulk} \oplus V_{\rm boundary}` (the bulk
  :term:`angular flux` :math:`\oplus` a single boundary trace — inflow and
  outflow fold into one :class:`~orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux`),
  posed :math:`A\psi = \tfrac{1}{k}F\psi` or :math:`A\psi = q`. The
  fission gain :math:`F` **never enters** :math:`A` — it is applied on
  the right-hand side and divided by :math:`k` at the eigenvalue layer.
  The boundary reflection is **no longer re-applied inside the streaming
  sweep** (the deleted "keystone"); it is delivered as the off-diagonal
  :math:`-B` coupling, and the outer Krylov / SI loop drives the
  boundary-consistency residual to zero. See :ref:`bc-extraction`.

- **Every operator's** ``.apply`` **output boundary is a**
  :class:`~orpheus.transport.source_sinks.angular_boundary_source_sink.AngularBoundarySourceSink`
  (Issue #208), completing
  the boundary half of the operator-output "dimensional-sin" carve (the
  bulk half — ``.apply.bulk`` →
  :class:`~orpheus.transport.source_sinks.angular_source_sink.AngularSourceSink`).
  The governing principle: an operator output
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
  (Issue #208/#201). The flux
  states :class:`~orpheus.transport.fields.angular_flux.AngularFlux` /
  :class:`~orpheus.transport.fields.scalar_flux.ScalarFlux` /
  :class:`~orpheus.transport.fields.harmonic_moment_flux.HarmonicMomentFlux`
  / :class:`~orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux` are an
  **affine space** :math:`\mathbb{A}` over a distinct **difference vector space**
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
  it**. A carrier is a cell
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
  :math:`C^1_{\rm int}` (Issue #208): the 2-D wavefront sweep and matvec no
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
  **retired** (#222): the cochain now lives in the rolling
  front (``_MovingFrontier``) and the full-cochain oracle history
  (``_octant_face_cochain``). See :ref:`wavefront-flux-cochain` for the
  succession.

- **The 2-D Cartesian SI iterate lives in moment space** — the
  within-group source-iteration fixed point
  :math:`\psi_{k+1} = (L{+}C)^{-1}(S\psi_k + B\psi_k + q)` consumes
  :math:`\psi` only through its flux moments :math:`\phi_\ell^m =
  (M\psi)_\ell^m`, so the *persistent* iterate is held as the moment
  tensor :class:`~orpheus.transport.fields.harmonic_moment_flux.HarmonicMomentFlux`
  (:math:`N \to (L{+}1)(2L{+}1)`, measured **18.3×** shrink at
  :math:`N=110`, :math:`L=1`) rather than the full :term:`per-ordinate <ordinate>`
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
  (Issue #208): the
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

- The transport eigenvalue and fixed-source problems share **one**
  operator algebra: they differ only in what sits on the right of the
  within-group transport operator :math:`A = L + C - S - B`.

  .. (vv-status rationale) The within-group fixed-source form
     A ψ = q, A = L+C-S-B. Verified by the operator-algebra assembly
     (build_within_group_system) and the fixed-source solver suites.
  .. vv-status: operator-fixed-source documented

  .. math::
     :label: operator-fixed-source

     A\,\psi \;=\; q
     \qquad \text{(fixed source)}

  .. (vv-status rationale) The k-eigenvalue form A ψ = (1/k) F ψ,
     A = L+C-S-B, with F the right-hand-side fission gain. Verified by
     the eigenvalue engines (power iteration and K = A⁻¹F) against the
     closed-form k∞ oracle.
  .. vv-status: operator-eigenvalue documented

  .. math::
     :label: operator-eigenvalue

     A\,\psi \;=\; \tfrac{1}{k}\,F\,\psi
     \qquad \text{(eigenvalue)}

  Both are built from operator addition, subtraction, scalar
  multiplication, and composition (``+``, ``-``, ``*``, ``@``) acting on
  :class:`~orpheus.numerics.operator.LinearOperator` instances; the
  fission gain :math:`F` is applied on the right and never enters
  :math:`A`.

- **The eigenvalue problem is layered into four tiers — leaves,
  posing, resolvent, algorithm**. The canonical **standard form** is the generalized
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
  :math:`A^{\dagger} = G^{-1} A^{\mathsf T} G`, NOT the Euclidean transpose. For the SN composite the
  metric :math:`G` is **block-diagonal** on :math:`V_{\rm
  bulk}\oplus V_{\rm trace}`: bulk :math:`V_{\rm cell}\,w_n` (phase-space
  measure) :math:`\oplus` trace :math:`|\Omega\cdot\hat n_f|\,w_n`
  (partial-current surface measure, pseudo-inverted on the singular
  tangential ordinates). The carrier is
  :class:`~orpheus.numerics.spaces.full_field_space.FullFieldSpace`;
  ``L``/``C``/``S``/``F``/``B`` all carry it so the
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
  :ref:`g-adjoint`.

- The :class:`~orpheus.numerics.operator.LinearOperator` Protocol
  carries one mandatory method (``apply``) and, per optional axis
  (inverse, adjoint, and **assembly**,
  :ref:`operator-algebra-assembly-axis`), a **three-layer structural
  surface** (#226): a runtime **predicate**
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
  ``phi.shape == (ng, nx, ny)`` for :term:`scalar flux`.  The canonical
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
  (GH #280/#282). The within-group
  augmented S\ :sub:`N` problem is posed as
  :math:`\bigl[\begin{smallmatrix} A_{AA} & A_{AB} \\ A_{BA} & A_{BB}
  \end{smallmatrix}\bigr]\bigl[\begin{smallmatrix} \psi_A \\ \psi_B
  \end{smallmatrix}\bigr] = \bigl[\begin{smallmatrix} q_A \\ q_B
  \end{smallmatrix}\bigr]` over **System A** (the transport
  :class:`~orpheus.transport.full_field.FullField`) and **System B** (the
  ψ½ radial-characteristic ray, its own
  :class:`~orpheus.transport.radial_characteristic_field.RadialCharacteristicField`
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

.. (vv-status rationale) The WDD cell-balance identity S·ψ̄ = source +
   streaming numerator, S = S_stream + σ_t V — the source of the
   σ-affine split. Definitional; the discrete sweep it feeds is verified
   downstream (dd-slab-scalar / dd-curvilinear-scalar).
.. vv-status: streaming-action-cell-balance documented

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
:ref:`sn-direct-seed-solve`.) The lesson
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
(:meth:`~orpheus.sn.operators.streaming.StreamingCollisionOperator.apply` still computes
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
:class:`~orpheus.sn.operators.streaming.StreamingCollisionOperator` (:math:`= L+C`), and
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
from :math:`A^{-1}` and :math:`B^{-1}`.

.. _smw-low-rank-exception:

The lone systematic exception is **low-rank structure**:
Sherman–Morrison–Woodbury rebuilds :math:`(A + U V^{\mathsf T})^{-1}`
from :math:`A^{-1}` at the price of one dense solve of rank size — and
that exception must be scoped honestly **per block**, not read as a
claim about the whole algebra (Issue #299). The dense collision
diagonal :math:`C` and the streaming operator :math:`L` carry no
low-rank split, so SMW buys nothing for the bulk pair. But the
boundary operator :math:`B` **is** exactly that structure: rank-1 per
face for the isotropically re-entering white/albedo laws (one
re-entry mode fed by a cosine-weighted outflow average), and an
ordinate permutation — rank :math:`N/2` on a slab face, trace-sized
rather than bulk-sized — for specular reflection (see the
:ref:`boundary-law census <bc-law-layer>` and
:ref:`bc-rank-n-algebra`). CP already exploits precisely this: its
white-boundary re-entry closes in **closed form** as the rank-1
update :math:`P_\infty = P_{\rm cell} + P_{\rm out} \otimes
P_{\rm in} / (1 - P_{\rm inout})` (``orpheus/cp/solver.py``), and
borrowing the same Woodbury closure for SN's boundary cycle — which
source iteration currently *iterates* — is Issue #300.

What each separate inverse would mean physically
------------------------------------------------

The within-group transport balance is the operator equation
:math:`(L+C)\,\psi = q`, i.e.

.. (vv-status rationale) The within-group transport balance
   Ω·∇ψ + Σ_t ψ = q. Governing equation (definitional).
.. vv-status: apply-solve-within-group-balance documented

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
global coupling enters. The :term:`diamond-difference <diamond difference>` cell update solves the
balance :eq:`streaming-action-cell-balance` for the cell-average flux:

.. (vv-status rationale) The per-cell WDD resolvent
   ψ̄ = (Q V + inflow)/S — the single-cell shadow of the non-distributing
   inverse. Definitional restatement; the sweep it expresses is verified
   downstream via dd-slab-scalar.
.. vv-status: apply-solve-cell-resolvent documented

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

  .. (vv-status rationale) The single-cell non-separability
     1/(S_stream + σ_t V) ≠ 1/S_stream + 1/(σ_t V). Mathematical
     identity, matching the sentineled apply-distributes /
     solve-does-not-distribute siblings.
  .. vv-status: apply-solve-denominator-inequality documented

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

.. (vv-status rationale) The operator-splitting (Neumann) series for
   (L+C)⁻¹ around the collision diagonal C. Mathematical identity
   (Lewis & Miller §3.2).
.. vv-status: apply-solve-neumann-series documented

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

.. (vv-status rationale) The term-by-term Neumann expansion
   C⁻¹ − C⁻¹LC⁻¹ + …. Mathematical identity (the expanded form of
   apply-solve-neumann-series).
.. vv-status: apply-solve-neumann-expansion documented

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

.. (vv-status rationale) The parallel / harmonic identity
   (L⁻¹ + C⁻¹)⁻¹ = L(L+C)⁻¹C. Mathematical identity (still not the
   coupled inverse).
.. vv-status: apply-solve-parallel-identity documented

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

.. (vv-status rationale) The source-iteration / Peierls
   collision-number series ψ = Σ_k [(L+C)⁻¹S]^k (L+C)⁻¹ q. Mathematical
   identity; its ρ = c convergence (green-scattering-ratio-bound) and the
   SI fixed point are exercised downstream by the source-iteration suites.
.. vv-status: apply-solve-source-iteration-series documented

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
\Sigma_{s,g}/\Sigma_{t,g} = c` (the :term:`scattering ratio`;
:ref:`affine-typed-field-algebra` documents the matching contraction
:math:`M = (L+C)^{-1}(S+B)` carried as a typed diagnostic on the SI
iterate). The sweep :math:`(L+C)^{-1}` being a *single bundled* inverse —
not :math:`L^{-1} + C^{-1}` — is exactly the point: it is the WDD
forward-substitution on :eq:`apply-solve-cell-resolvent`, dividing by the
summed denominator cell-by-cell in inflow-to-outflow order. See Lewis &
Miller, *Computational Methods of Neutron Transport* (:cite:`LewisMiller1984`,
§3.2 for the sweep as the discrete-ordinates resolvent and §4 for the
source-iteration / Neumann scattering series), and Adams & Larsen 2002
(:cite:`AdamsLarsen2002`, §II for the spectral radius :math:`\rho = c`).

Why this is the right architecture, not a limitation
----------------------------------------------------

Invertibility is a property of the **sum**, not of the parts. That is
exactly why :math:`L + C` is packaged as one
:class:`~orpheus.sn.operators.streaming.StreamingCollisionOperator`: the
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
  :class:`~orpheus.sn.operators.streaming.StreamingCollisionOperator` reports
  :attr:`~orpheus.numerics.operator.LinearOperator.is_invertible`
  ``= True`` and carries a direct-sweep ``solve``; the leaves do not
  (streaming has no ``solve`` at all; collision's ``solve`` is the
  *local* :math:`q/\sigma`, which is the
  :math:`C^{-1}` of a *different* problem, never the coupled inverse).
  The :class:`~orpheus.numerics.operator.OperatorSum` deliberately
  **does not propagate** ``solve`` (:ref:`composition-algebra`); the
  :class:`~orpheus.sn.operators.streaming.StreamingCollisionOperator` *adds it back* via the
  SN-specific algebraic identity "WDD sweep :math:`\approx (L+C)^{-1}`"
  (:cite:`LewisMiller1984` §3.2). The composite owns ``apply``, ``solve``, and
  ``apply_transpose`` as three actions of **one** operator on a single
  shared
  :class:`~orpheus.sn.loss_representation.LossRepresentation` (L21 —
  "matvec ≡ sweep"); :meth:`StreamingCollisionOperator.apply` and
  :meth:`StreamingCollisionOperator.solve` single-source :math:`\sigma` from the
  collision diagonal, so they cannot disagree on which loss they invert.

The :ref:`three-layer operator surface <capability-set-semantics>` is
what makes this architecture *enforced* rather than merely *intended*: a
downstream Krylov consumer that asks for the inverse of a bare
:class:`~orpheus.sn.operators.streaming.StreamingOperator` cannot even
spell ``L.inverse()`` (the streaming leaf declares no such method — a
*static* error), and a generic sum that has not been promoted to
:class:`~orpheus.sn.operators.streaming.StreamingCollisionOperator` returns the
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

.. (V&V scope note) The inverse-as-operator keystone (#226 Phase 2):
   applying an operator's inverse OBJECT equals invoking its native
   realization verb ``solve``. A foundation software invariant (no
   eigenvalue / flux claim); the label is wired to the bit-identity
   gate ``tests/sn/operators/test_inverse_operator_equivalence.py``
   (``(L+C).inverse().apply(b) == (L+C).solve(b)`` for the sweep-invertible
   loss operator, plus the seed-drop and returned-surface-type checks).

.. math::
   :label: inverse-as-operator

   A^{-1} b \;=\; \texttt{A.inverse().apply(b)} \;=\; \texttt{A.solve(b)}

For the sweep-invertible loss operator :math:`(L+C)` the native ``solve``
IS the WDD sweep, so this keystone reads: applying the inverse object runs
the same sweep the operator's own ``solve`` runs — no separate inverse
machinery.

**Deleted** (the algebra-closed and driver-realized kinds):

- :class:`~orpheus.numerics.operator.OperatorSum` — a generic sum's
  inverse action is driver-realized by the
  :class:`~orpheus.numerics.green_operator.GreenOperator` (the
  preconditioned splitting), not a substrate verb. (The
  sweep-invertible ``(L+C)``
  :class:`~orpheus.sn.operators.streaming.StreamingCollisionOperator` subclass
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
       structure — which the boundary block :math:`B` HAS, unlike
       the bulk :math:`C`, :math:`L`;
       :ref:`the scoped statement <smw-low-rank-exception>`).
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

When a composed within-group operator meets the Krylov accelerator's
:mod:`scipy.sparse.linalg` interface, the ORPHEUS↔scipy boundary is
crossed at a **single** site internal to
:class:`~orpheus.numerics.iteration.KrylovAcceleration`: one ravel-aware
adapter ``_as_scipy_linop(carrier_matvec, template, n)`` lifts the flat
scipy vector to the typed carrier (via ``from_flat``), runs the
carrier-space matvec, and ravels the result back — so both the system
matvec and the preconditioner route through one source of truth (#257
S7; the per-operator ``build_transport_linear_operator`` scipy wrappers
it replaced are retired). The system matvec is the named carrier-space
closure ``loss_minus_gains`` — the invertible resolvent minus the lagged
coupling gains, :math:`(L{+}C)\,\psi - \sum_{\rm gains} g\,\psi`, i.e.
the honest full within-group operator :math:`A = L + C - S - B` applied
(the gains are the scattering :math:`S` and the boundary reflection
:math:`B`) — which reads like the operator rather than ravel plumbing.
Because it is expressed purely through the operator algebra, it is the
discretisation-agnostic form a unified cross-solver iteration driver
(SN / MoC / CP / diffusion, Issue #14) consumes without knowing which
transport discretisation produced the triple. A future consumer needing
a standalone operator→scipy adapter (e.g. a DSA preconditioner built as
an operator, Issue #2) should **generalise** ``_as_scipy_linop`` to
accept an ``op.apply`` callable, **not** resurrect the retired flat-only
``as_scipy_linop`` twin.


.. _intrinsic-operator-types:

The intrinsic operator types
============================

The maps of the algebra fall into a three-way partition by what their
action *returns* — the grand report §5.6 **suffix law**. An
**Operator** carries a field to a field: the local multiplication
:math:`C = M[\sigma_t]` (the collision diagonal, §5.7) and the
streaming leaf :math:`L` are its transport instances. A **Kernel** is
the nonlocal refinement — an Operator whose output at a point
*integrates* the input across an axis — realized by the two Boltzmann
emission kernels, scattering :math:`S = R\circ\Lambda\circ M` and
fission :math:`F = |\chi\rangle\langle\nu\Sigma_f|`. A **Functional**
is the disjoint sibling that maps a field to a **scalar** — the
reaction-rate integrals behind the :math:`k`-eigenvalue. The
discriminator between the local Operator and the nonlocal Kernel is
**locality**: multiplication is the diagonal sub-algebra, whose output
at a phase-space point reads the input only there, and the kernels are
everything off-diagonal. The three type sections below develop this
partition in order — the diagonal / multiplication Operator, then the
Functional, then the two integral Kernels — and the full codomain
partition with its type-system table is set out in
:ref:`functional-category`.


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
:class:`~orpheus.sn.operators.streaming.StreamingCollisionOperator` lives **one-directionally**
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
   :class:`~orpheus.sn.operators.streaming.StreamingCollisionOperator` has its own
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
   diffusion arm — see :doc:`/theory/methods/diffusion_1d`), but poses its
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
   :doc:`/theory/methods/sn/index`.

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
(:doc:`/theory/foundations/frame`).

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
   kernel (the kernel is the redistribution, not the :term:`quadrature`
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

.. (vv-status rationale) The carrier-grid composite
   S_aniso = (1/W)(R∘Λ∘M). Representational (named-field-typing)
   identity, matching the sentineled scattering-carrier-grid /
   carrier-grid-* siblings; the bit-identical composed kernel is pinned
   by test_scattering_kernel_crosscheck.
.. vv-status: scattering-aniso-composite documented

.. math::
   :label: scattering-aniso-composite

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

.. (vv-status rationale) The Liskov typing identity HarmonicFrame IS-A
   GalerkinFrame (the angular SH projection is the canonical
   pure-Galerkin frame). Structural / representational identity; the
   from_galerkin faces are bit-identical to the generic numerics frame's.
.. vv-status: harmonic-frame-is-galerkin documented

.. math::
   :label: harmonic-frame-is-galerkin

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
   0-ULP equivalence). Deferred follow-up: #260
   (:class:`~orpheus.numerics.operator.SumOfTensorProductsOperator`
   un-orphaning).


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


Boundary conditions as operators
================================

Every *realised* boundary condition IS a Wave-0
:class:`~orpheus.numerics.operator.LinearOperator`, composable with the
streaming / collision / scattering / fission operators through the same
algebra dunders — which is why the boundary law :math:`B` enters
:math:`A = L + C - S - B` as a first-class sibling rather than being
folded into streaming. A boundary *law* is a pure descriptor with **no**
``apply``; only its *realisation* through
:class:`~orpheus.sn.boundary.realizer.SNBoundaryRealizer` is a Wave-0
operator.

The full treatment — the §15.2 :math:`G_\alpha` geometric realisation
primitives (permutation, mask, average, wrap, source; the law→primitive
map), the rank-N Marshak / partial-current **descriptor-tree algebra**
(``LawSum`` / ``LawScaled`` realised by
:func:`~orpheus.geometry.boundary.realize_recursively`), and the
load-bearing **descriptor-tree vs operator-tree** type separation
("mixing the two algebras is a type error") — is documented on the
boundary-condition page: see the :ref:`extraction narrative
<bc-extraction>`, the :ref:`law→primitive realisation map
<bc-tensor-primitives>`, and the :ref:`recursive realiser
<bc-realize-recursively>` in
:doc:`/theory/foundations/boundary_conditions`.

Tensor-network decomposition
============================

The factored / tensor-network shape decomposition of the S\ :sub:`N`
operators — *which* algebraic shape (tensor product, sum-of-tensor-products,
or an irreducible sequential walk) each operator leaf takes, the MA-Q1
admissibility condition, and why streaming resists a clean factorisation —
is developed in :doc:`/theory/foundations/operator_tensor_network`.


The affine-typed field algebra
==============================

The affine-typed S\ :sub:`N` **field** algebra — the recognition that
flux **states** form an affine space :math:`\mathbb{A}` (points, no
origin) over a difference vector space :math:`V`, with the field-role
grid completed by a **displacement** column that carries the convergence
diagnostics a state cannot — is developed in full on its own page,
:doc:`/theory/foundations/field_algebra`. That page types the **fields**
(state / displacement / residual / source) that the operators of this
page act on; it also wires the typed equation residual
:math:`r = (L + C - S - B)\,\psi - q` to its ``from_balance`` consumer.

The composite metric adjoint
============================

The **metric-correct Hilbert adjoint** ``op.H`` of the composed loss
operator :math:`A = L + C - S - B` — the **G-adjoint**
:math:`A^{\dagger} = G^{-1} A^{\mathsf T} G` over the block-diagonal
:class:`~orpheus.numerics.spaces.full_field_space.FullFieldSpace` metric,
with the singular trace block inverted by a Moore–Penrose pseudo-inverse
— is developed on its own page, :doc:`/theory/foundations/operator_adjoint`.
That page derives it from the reciprocity identity :math:`\langle A\psi,
\varphi\rangle_G = \langle \psi, A^{\dagger}\varphi\rangle_G`, shows why
the metric applies once at the operator level, and pins reciprocity to
round-off against a structurally-independent dense-probe oracle. It is the
adjoint face of the operator algebra on this page — distinct from the
frame's Petrov–Galerkin test-space adjoint in
:doc:`/theory/foundations/frame`.


The interior face-flux cochain
==============================

The **interior face-flux cochain** :math:`C^1_{\rm int}` — the
sweep-internal record of angular flux on interior faces, its biproduct
decomposition :math:`C^1 = C^1_{\rm int} \oplus C^1_\partial` and trace
algebra, and the succession note on why the typed ``WavefrontFlux``
carrier retired (the concept survives in its two native realizations) —
is developed in full on its own page,
:doc:`/theory/foundations/wavefront_cochain`.

The inverse family
==================

How the composed operator's inverse is realized and materialized — the
driver-applied sweep, the Green preconditioned inverse, the dense
materialising inverse, and the sparse assembly axis — is developed in
:doc:`/theory/foundations/operator_inverse_family`.

The coupled block operator
==========================

The **coupled block operator** — the 2×2 block system in which the
curvilinear starting-direction flux :math:`\psi_{1/2}` (the ψ½ ray)
becomes a first-class **System B** alongside the angular-flux System A,
its four named blocks (:math:`A_{AA}` the within-group loss composite,
plus the seed / emission / march couplings), the N-general block
machinery, and the structure-keyed block solve — is developed in full
on its own page, :doc:`/theory/foundations/coupled_block_operator`.

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

The two halves are directional selectors on the single unified
:class:`~orpheus.numerics.space.FunctionSpace` subclass
:class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace`
(the per-face ``InflowTraceSpace`` / ``OutflowTraceSpace`` classes were
consolidated into inflow/outflow selectors on it — the
:meth:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.inflow_indices_for_face`
/ :meth:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.outflow_indices_for_face`
methods), which carry a **per-face directional mask**:

.. math::
   :label: per-face-inflow-mask

   \mathrm{inflow\_mask}[f, n]
   \;=\;
   \bigl(\Omega_n \cdot \hat n_f < -\epsilon\bigr).

.. vv-status: per-face-inflow-mask documented

The mask has shape ``(n_faces, n_ordinates)`` boolean; the
tangential band :math:`|\Omega_n \cdot \hat n_f| \leq \epsilon`
(``TANGENTIAL_EPS = 4 * np.finfo(np.float64).eps``, i.e.
:math:`\approx 8.9\times 10^{-16}`) is in neither half — so **"not
inflow" is never "outflow"**, and the two index sets are disjoint but
NOT exhaustive.

Three consumers read the mask today:

* The SN realizer consumes **both** selectors: ``inflow_indices_for_face(face)``
  gives a realized law's codomain :math:`\Gamma_-`, and
  ``outflow_indices_for_face(face)`` its **domain** :math:`\Gamma_+` —
  restricted through the
  :class:`~orpheus.numerics.operator.TraceRestrictionOperator` the same
  table caches (campaign phase B3.2; see
  :ref:`bc-domain-narrowing`).
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

.. (vv-status rationale) The resolvent standard form
   K_pm ≡ A_loss⁻¹ M with K_pm ψ = (1/λ)ψ. Definitional / posing
   identity, matching the sentineled eigen-standard-form; the
   k = λ_max(A⁻¹F) claim is anchored by the homogeneous closed-form
   algebra-of-record (kinf_and_spectrum_homogeneous).
.. vv-status: eigen-resolvent documented

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
       every posing (step O.2b R5, :ref:`g-adjoint`).
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
:class:`~orpheus.numerics.iteration.KEigenvalue` (built from the
``(A, S, F)`` triple) is the
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

.. (vv-status rationale) The k-row posing (L+C−S−B)ψ = (1/k)Fψ ⟺
   [(L+C−S−B)⁻¹F]ψ = kψ — eigen-standard-form specialized to the LIVE
   k-row. Definitional / posing identity; k_eff is anchored by the
   homogeneous closed-form oracle (kinf_and_spectrum_homogeneous).
.. vv-status: eigen-k-posing documented

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

.. (vv-status rationale) The α-eigenvalue derivation (the e^{αt} ansatz
   → (L+C−S−F−B)ψ = −α T ψ). Definitional derivation of a documented
   future seam — the α-row is Not built (only the k-row exists;
   unify-after-two).
.. vv-status: eigen-alpha-derivation documented

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

**The adjoint row (LIVE — #276 A4/A5).** The adjoint eigenproblem
:math:`A_{\rm loss}^{\dagger}\,\psi^{\dagger} = \lambda\,M^{\dagger}\,
\psi^{\dagger}` is **just another posing row** whose role-operators are
the daggers of the forward leaves — and it now RUNS in production:
``KEigenvalue((L+C).H, (S+B).H, F.H)`` through the unchanged
:func:`~orpheus.numerics.eigenvalue.power_iteration`
(:func:`~orpheus.sn.solver.solve_sn_adjoint`; the full chapter is
:ref:`sn-adjoint`). The dagger is *free* from the
dagger-biproduct category already documented on this page (the ``.H``
adjoint propagates through ``+`` / ``&`` / ``@`` — see
:eq:`tensor-product-adjoint-distributivity`) and, crucially, from the
**composite metric-correct G-adjoint** that step O.2b R5 made concrete:
``op.H`` is the metric-correct :math:`A^{\dagger} = G^{-1} A^{\mathsf T}
G` over the :math:`V_{\rm bulk}\oplus V_{\rm trace}` carrier
(:ref:`g-adjoint`), dense-oracle-verified to round-off.
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
:doc:`/theory/methods/diffusion_1d`), yet it still belongs to the monolithic
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
:ref:`g-adjoint`): a block-diagonal Hilbert metric with a
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
(:ref:`g-adjoint`). Since P4.5 W-D the bulk leaves
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
   * - 2026-07-12
     - **The curvilinear ψ½ ray is System B of a 2×2 coupled block
       operator** — the augmented S\ :sub:`N` within-group problem is
       posed as a 2×2 block system over **System A** (the transport
       :class:`~orpheus.transport.full_field.FullField`) and **System B**
       (the ψ½ radial-characteristic ray). The four blocks — the seed
       :math:`A_{AB}`, the emission :math:`A_{BA}`, the radial march
       :math:`A_{BB}` (its direct Carlson solve IS the inverse) — are
       named operators assembled by the N-general
       :class:`~orpheus.numerics.coupled_system.CoupledOperator`
       machinery; :func:`~orpheus.sn.coupled_system.build_within_group_system`
       is the one production spelling, and System B's presence is
       **structural** (it exists iff the mesh carries a ray). See
       :ref:`coupled-block-operator`.
     - #280 / #282
     - ``main`` (``6732778a``)
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
   * - 2026-06-07
     - **The 2-D Cartesian SI iterate lives in moment space** — the
       within-group source-iteration fixed point consumes :math:`\psi`
       only through its flux moments, so the *persistent* iterate is held
       as the moment tensor
       :class:`~orpheus.transport.fields.harmonic_moment_flux.HarmonicMomentFlux`
       rather than the full per-ordinate
       :class:`~orpheus.transport.fields.angular_flux.AngularFlux`; the
       source stays bit-identical (shared :math:`R\,\Lambda`
       reconstruction), only the SI convergence test moves to the moment
       :math:`L^2` (principled-equivalence). 2-D Cartesian, interior-bulk
       only. See :ref:`sn-angular-windowing`.
     - #205
     - ``main`` (Phase 5a, ``93807aa`` / ``b97d4f9`` / ``13ca001``)
   * - 2026-06-05
     - **The eigenvalue problem is layered into four tiers** (leaves /
       posing / resolvent / algorithm): the generalized eigenproblem
       :math:`A_{\rm loss}\,\psi = \lambda\,M\,\psi` (k-row
       :math:`A_{\rm loss} = L+C-S-B`, :math:`M = F`) has its power-method
       realization as the dominant eigenpair of the resolvent
       :math:`A_{\rm loss}^{-1} M`;
       :func:`~orpheus.numerics.eigenvalue.power_iteration` is the
       canonical Layer-4 algorithm and
       :class:`~orpheus.numerics.iteration.KEigenvalue` delegates its loop
       to it (one loop, Cardinal Rule 2). See :ref:`eigenvalue-posing`.
     - —
     - ``main`` (``650032e`` / ``7603c8e``)
   * - 2026-06-05
     - **The Hilbert adjoint** ``op.H`` **is the metric-correct
       G-adjoint** :math:`A^{\dagger} = G^{-1} A^{\mathsf T} G` over the
       block-diagonal
       :class:`~orpheus.numerics.spaces.full_field_space.FullFieldSpace`
       metric (:math:`V_{\rm bulk} \oplus V_{\rm trace}`; singular trace
       block pseudo-inverted) — NOT the Euclidean transpose; applied
       **once at the operator level** and raising
       :class:`~orpheus.numerics.operator.MissingAdjoint` eagerly rather
       than silently going Euclidean. See :ref:`g-adjoint`.
     - —
     - ``main`` (Wave O step O.2b R5, ``5c06196``)
   * - 2026-06-04
     - **2-D Cartesian eigenvalue problems solve via BOTH inner solvers**
       — the source-iteration inner is the geometry-agnostic structural
       twin of the Krylov inner: identical composite RHS, identical loss
       decomposition (the invertible resolvent :math:`L + C` plus the two
       lagged coupling gains — scattering :math:`S`, boundary reflection
       :math:`B` — on the variadic driver), differing **only** in the
       iteration driver. See :ref:`bc-extraction-2d-si-krylov-twin`.
     - #208
     - ``main`` (Wave O 2-D SI Phase A)
   * - 2026-06-04
     - **The interior cell-face angular fluxes are a 1-cochain**
       :math:`C^1_{\rm int}` — with the boundary trace
       (:math:`C^1_\partial`) it biproduct-decomposes the full face
       cochain :math:`C^1 = C^1_{\rm int} \oplus C^1_\partial`; the
       sweep's seed/absorb are the typed trace operators :math:`\iota_*`
       / :math:`\iota^*` (``absorption = identity`` is the biproduct law
       :math:`\iota^* \circ \iota_* = \mathrm{id}`). The ``WavefrontFlux``
       carrier retired at S6.4(f). See :ref:`wavefront-flux-cochain`.
     - #208 / #222
     - ``main`` (Wave O #205 Phase 5)
   * - 2026-06-03
     - **The boundary law** :math:`B` **becomes a first-class sibling
       operator** and every operator's ``.apply`` output is typed as a
       *source/sink* (the bulk → :class:`~orpheus.transport.source_sinks.angular_source_sink.AngularSourceSink`,
       the boundary → :class:`~orpheus.transport.source_sinks.angular_boundary_source_sink.AngularBoundarySourceSink`),
       completing the operator-output role typing
       (:ref:`bc-extraction-operator-output-typing`).
     - #208
     - ``main`` (Wave O steps O.4a.2 / B.5.2)

