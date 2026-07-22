.. _affine-field-algebra:

==============================
The Affine-Typed Field Algebra
==============================

.. contents:: Contents
   :local:
   :depth: 2


.. Machine header — the ``nexus-meta`` schema for this page (PROVISIONAL).
.. This page was extracted from ``operator_algebra.rst`` (#231 Phase 3); its
.. content is verbatim except for the split-time symbol reservation noted
.. below. The schema is provisional pending a full re-audit of the corpus.

.. dropdown:: Machine header — ``nexus-meta`` schema (PROVISIONAL)
   :color: muted

   .. code-block:: yaml

      module: transport
      concept: field_algebra
      role: "the affine-typed SN flux field algebra — flux states as an affine torsor over the displacement space, and the state / displacement / residual / source role grid"
      depends_on: [operator_algebra]
      related: [operator_adjoint, boundary_conditions]
      status: "extracted from operator_algebra.rst; content verbatim, provisional header"


This page develops the **affine-typed field algebra** of the S\ :sub:`N`
transport solve — the recognition that flux **states** do not form a
vector space but an **affine space** :math:`\mathbb{A}` over a difference
(displacement) vector space :math:`V`, and that the field-role grid is
completed by a fourth **displacement** column carrying the convergence
diagnostics a state structurally cannot. It is the field-algebra companion
to the operator algebra developed in
:doc:`/theory/foundations/operator_algebra`: that page types the
**operators** — the within-group loss composite :math:`A = L + C - S - B`
and its invertible sub-composite :math:`L + C`, whose inverse
:math:`(L+C)^{-1}` is the transport :term:`sweep` — while this page types the
**fields** those operators act on.

.. note::

   **Two symbols that look alike are kept distinct throughout this page.**
   :math:`\mathbb{A}` (blackboard-bold) is the **affine space of flux
   states** — points with no origin, over the displacement space
   :math:`V`. :math:`A = L + C - S - B` (italic) is the **within-group
   loss operator**. They are different objects: the operator acts *on*
   points of the affine space. The **sweep** is the inverse of the
   invertible sub-composite, :math:`(L+C)^{-1}` — the *inner kernel* of
   the full within-group solve :math:`A^{-1}`, never :math:`A^{-1}`
   itself.

.. _affine-typed-field-algebra:

The affine-typed SN field algebra — state, displacement, residual
=================================================================

Wave O step **O.2** (`Issue #208
<https://github.com/deOliveira-R/ORPHEUS/issues/208>`_, which subsumes
`Issue #201 <https://github.com/deOliveira-R/ORPHEUS/issues/201>`_,
2026-06-08) closes the operator-algebra-completeness wave by making the
SN **flux algebra** as type-honest as the operator algebra itself.
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
     :math:`A = L + C - S - B`. Flux and rate-density are joined by
     :math:`A` (:math:`A`: apply, flux → rate-density;
     :math:`A^{-1}`: the full within-group solve, rate-density → flux).
     The invertible sub-composite :math:`L + C` has the **sweep**
     :math:`(L+C)^{-1}` as its inverse — the inner kernel of
     :math:`A^{-1}`, NOT :math:`A^{-1}` itself. State and displacement
     live in the **flux** universe; source and residual live in the
     **rate-density** universe.
   - **Flux states are an affine space** :math:`\mathbb{A}` **over a
     difference vector space** :math:`V`. The increment
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
operator :math:`A = L + C - S - B`:

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

The map between them is dimensional: applying :math:`A` (the full
within-group loss — streaming :math:`\hat\Omega\cdot\nabla` + collision
:math:`\Sigma_t` − scattering − boundary, each carrying a
:math:`1/\mathrm{cm}` factor) sends a flux to a rate-density (``A.apply`` →
:class:`AngularSourceSink`, documented at
:ref:`bc-extraction-operator-output-typing`); inverting it
(``A.solve`` :math:`= A^{-1}`, the full within-group solve) sends a
rate-density back to a flux (``.solve`` → :class:`AngularFlux`). The
invertible sub-composite :math:`L + C` has the sweep :math:`(L+C)^{-1}`
as its inverse — the inner kernel that ``A.solve`` drives to its fixed
point, never :math:`A^{-1}` itself. The source side (the B.5.2
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
space :math:`\mathbb{A}` has points but no distinguished origin; differences of
points live in an associated **difference vector space** :math:`V`
(the *torsor* :math:`\mathbb{A}` over :math:`V`). The torsor axioms ARE the
field algebra (the SymPy-free derivation lives in
:class:`~orpheus.transport.fields._flux_role.FluxRole`):

.. math::
   :label: affine-torsor-algebra

   \psi_2 \ominus \psi_1 &\;\to\; \Delta\psi \in V
       \qquad&&\text{(the iterate increment; the ONLY mint of a displacement)}\\
   \psi \oplus \Delta\psi &\;\to\; \psi' \in \mathbb{A}
       \qquad&&\text{(the torsor action } \mathbb{A}\times V\to \mathbb{A}\text{; the update step)}\\
   \psi_1 \oplus \psi_2 &\;\to\; \bot
       \qquad&&\text{(no origin — } \textstyle\sum\lambda_i = 2 \text{ lands off } \mathbb{A})\\
   \textstyle\sum_i \lambda_i \psi_i,\ \textstyle\sum_i\lambda_i=1
       &\;\to\; \psi \in \mathbb{A}
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
  the affine subspace :math:`\mathbb{A}`. The :class:`TypeError` names the legal
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
adding two of them is meaningless (it lands off :math:`\mathbb{A}`). The same
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
     - flux states are points in :math:`\mathbb{A}`; :math:`\Delta\psi` is a
       vector in :math:`V`; the three ops are the torsor axioms
       :eq:`affine-torsor-algebra`
     - ``state + state`` becomes UNREPRESENTABLE by construction; the
       #201 gate is a type consequence, not a runtime check
   * - **Banach fixed-point / contraction**
     - :math:`\Delta\psi^{(i+1)} = M\,\Delta\psi^{(i)}`,
       :math:`M = (L+C)^{-1}(S+B)` — an exact linear recurrence on the
       increments
     - the displacement is the natural carrier of the contraction
       factor :math:`\rho`, the a-posteriori bound, the Aitken
       extrapolation — data a state has nowhere to put
   * - **Krylov / residual dual**
     - :math:`\Delta\psi = (L+C)^{-1}\,\tilde r` — the increment and the
       residual are the SAME defect, mapped between the two universes by
       the sweep :math:`(L+C)^{-1}` (the inner kernel of :math:`A^{-1}`)
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
:math:`M = (L+C)^{-1}(S+B)`, a contraction with spectral radius
:math:`\rho(M) \le \max_g \Sigma_{s,g}/\Sigma_{t,g} = c`
(:cite:`AdamsLarsen2002`).

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
true error by a factor :math:`1/(1-\rho)`. As the :term:`scattering ratio`
:math:`c \to 1` (optically thick, near-pure-scatter), :math:`\rho \to
1` and the understatement blows up: at :math:`c = 0.99`,
:math:`\rho\approx 0.99` and the true error is :math:`\sim 100\times`
the increment — so a solve that "converges" at
:math:`\lVert\Delta\psi\rVert < \text{tol}` is actually
:math:`\sim 100\cdot\text{tol}` from the solution. This is the
canonical source-iteration stall-masking-as-convergence trap
(:cite:`AdamsLarsen2002`).
:meth:`~orpheus.transport.displacements._displacement.Displacement.true_error_estimate`
surfaces it (it raises if :math:`\rho \notin [0,1)` — a
non-contracting iteration has no finite geometric-tail estimate). The
displacement is the natural home for this because the estimate needs
both the increment AND :math:`\rho` (which needs the *previous*
increment) — data a flux state cannot hold.

**The convergence map.**
:meth:`~orpheus.transport.displacements._displacement.Displacement.where_largest`
returns the :math:`k` index tuples with the largest
:math:`|\Delta\psi|` — the per-cell / per-group / :term:`per-ordinate <ordinate>`
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
:doc:`/theory/methods/diffusion_1d`); :eq:`affine-typed-residual-eq` is that
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
