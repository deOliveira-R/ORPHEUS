.. _sn-solver-operator-algebra-coordinator:

SNSolver as an operator-algebra coordinator
============================================

This chapter is where the book's machinery converges: the operators
the preceding chapters built — the swept :math:`(L+C)`, the scattering
:math:`S`, the boundary gain :math:`B`, the fission :math:`F` —
coordinated into the production eigenvalue solve, and, from the
converged flux, the frame projections (homogenisation, condensation)
that hand a coarse problem back to the same solver.

⛔ **Since 2026-09-13,** :class:`~orpheus.sn.solver.SNSolver` **caches no
operator at all.**  It held three until that day, and lost all three to
the consumers campaign's step 2 in the order the campaign's own units
landed:

.. list-table::
   :header-rows: 1
   :widths: 26 20 54

   * - retired slot
     - went where
     - what reads it now
   * - ``SNSolver.fission_op``
     - the Problem hub (unit C2)
     - :attr:`SNProblem.fission <orpheus.sn.problem.SNProblem.fission>`
       — the *composite*
       :class:`~orpheus.transport.operators.fission.FissionOperator` on
       the full field.
       :meth:`~orpheus.sn.solver.SNSolver.compute_fission_source` applies
       its derived energy face ``problem.fission.isotropic_energy``,
       bit-identically to the retired mint
       (:ref:`sn-one-fission-per-problem`)
   * - ``SNSolver.scattering_op``
     - the Problem's posed record (unit C3b)
     - ``problem.system.factors.scattering`` — the
       :class:`~orpheus.transport.operators.scattering.ScatteringOperator`
       carrying the P0 in-scatter + the P\ :sub:`ℓ` Galerkin
       reconstruction (Wave D Issue 13), minted ONCE per Problem at the
       hub's retained Legendre order
   * - ``SNSolver.n2n_op``
     - the Problem's posed record (unit C3b)
     - ``problem.system.factors.n2n`` — the
       :class:`~orpheus.transport.operators.n2n.N2NOperator`.  It was a
       **passenger inside** ``scattering_op`` until CS4c step 3
       (2026-08-30), when the channel became first-class because its
       bundling — scattering-like or production-like — is
       context-dependent and must not be decided at the operator level
       (:ref:`n2n-reactions`, :ref:`sn-n2n-adjoint`).  Every
       ``(n,2n)`` verb on the solver — since #448 that is the group-rate
       accumulations alone; the ``_add_n2n_source`` delegator retired
       with the hand-built finalize source that was its only caller —
       routes through its energy binding's field

.. note:: **What the head of this chapter used to say, and why every
   version of it was a description of a seam.**

   Until 2026-09-13 it read *"the solver caches two operators … which
   read their cross sections through the hub's one*
   :class:`~orpheus.transport.mesh.material_xs_field.MaterialXSField`
   *rather than snapshotting values, and therefore survive a rebind
   untouched"*.  The trailing clause went first: there **is** no rebind
   any more — :math:`\sigma_t` is a datum of the Problem and a σ-variant
   is another Problem (:ref:`sn-sigma-is-a-problem-datum`).  Then the
   subject went too, because a solver-held copy of a Problem's leaf is
   the thing the campaign exists to remove.

   Read-through was never the property that mattered; it was the
   *mitigation* for holding a copy of something that was not the
   solver's.  The copies were also an injection **seam**: the builder
   took ``scattering_op=`` / ``n2n_op=`` keyword arguments so a caller
   could hand it a foreign :math:`S`, which meant the Problem's
   :math:`A` was not a function of the Problem's data.  ``[M]`` 34
   keyword injections across 14 files, of which **33** handed back the
   hub's own operator (the seam being paid for and not used) and **one**
   genuinely minted a foreign P0 :math:`S`.  Both keywords are
   **deleted**, and that one call site now poses on
   ``sn.system`` like everything else
   (:ref:`sn-the-problem-poses-its-pencil`).

The story that section closes began at CS4c step 4 (2026-08-30), when
the fission channel became *two bindings of one datum*
(:ref:`sn-fission-binding-adjoint`) and every consumer was re-pointed at
the one it actually feeds — the k-outer hands bare
:math:`(n_g, *\text{spatial})` scalar arrays, so it wants the energy
binding on the Problem's **bulk** space, while the eigen-:math:`M` posing
wants the composite.  Step 4 gave each consumer the right *binding* and
left each minting its own *object*; step 2 of the consumers campaign
made it one object with two faces.

.. note:: **The retained Legendre order is NOT a solver argument — it is
   the hub's datum.**

   Until 2026-09-12 :class:`~orpheus.sn.solver.SNSolver` took its own
   ``scattering_order`` keyword and clamped it; the adjoint entries and
   the fixed-source entry each had their own, and the three disagreed on
   a live solve.  Since the consumers campaign's step 1 (ruling R-cc9;
   GitHub #459) the order is clamped ONCE on
   :class:`~orpheus.sn.problem.SNProblem` at construction, and
   ``SNSolver.__init__`` reads :attr:`problem.scattering_order
   <orpheus.sn.problem.SNProblem.scattering_order>`.  The
   attribute ``SNSolver.scattering_order`` still exists and still means
   the same thing — it is now a **read**, not a second spelling.  The
   four entry
   points keep their ``scattering_order=`` keyword as sugar that
   forwards into the hub's constructor, so no user-facing call changed.
   The argument and the clamp are documented at
   :ref:`sn-hub-retained-order`.

The loss composite :math:`L+C` is the **Problem's**, and there is exactly
one of it: ``problem.system.factors.streaming_collision``, built through
the one spelling
:func:`~orpheus.sn.coupled_system.build_streaming_collision` inside the
one builder call, and it *is* the object every Strategy value inverts.

.. note:: **The in-code ruling this paragraph mirrors was rewritten on
   2026-09-13, and the rewrite is a strengthening rather than a
   reversal.**

   It used to read: :math:`L+C` is deliberately **not** cached on the
   solver, because a solver-held second copy would be a twin free to
   drift from the operand the sweep uses — and the twin hazard was not
   hypothetical, it had already happened (the retired
   ``SNSolver.L``/``.S``/``.F`` triple below).  A second clause made the
   drift concrete: a
   :class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`
   holds its coefficient as a **snapshot**, so a cached :math:`C` would
   go stale the moment :math:`\sigma_t` was rebound underneath it.

   Both halves are now discharged by a stronger move rather than by a
   prohibition.  The staleness half **dissolved**: there is no rebind,
   because :math:`\sigma_t` is a Problem datum
   (:ref:`sn-sigma-is-a-problem-datum`).  The twin half is **honoured**:
   the hub's record is the only copy, the solver holds none, and "do not
   cache it here" became "there is nothing here to cache".

.. note::

   Before 2026-07-28 the solver also exposed an ``(L, S, F)`` "operator
   triple" (``SNSolver.L`` / ``.S`` / ``.F``).  It was **retired**: the
   attributes had no production reader, and ``SNSolver.L`` was a
   misnomer — it held the *composite* :math:`L+C`, whereas :math:`L`
   throughout this book is the :math:`\sigma`-free streaming leaf
   (:ref:`the affine collision split <operator-algebra>`).  Consumers
   needing the composite read it off the Problem's factors; the
   free-function spelling
   :func:`~orpheus.sn.coupled_system.build_streaming_collision` is the
   builder's, not a consumer entry point.

Every leaf on the Problem's record is a
:class:`~orpheus.numerics.operator.LinearOperator`
in the Wave A operator-algebra sense: predicate-typed, composable
under :class:`~orpheus.numerics.operator.OperatorSum` and
:class:`~orpheus.numerics.operator.OperatorProduct`, and protocol-
conforming so the iteration primitives in
:mod:`orpheus.numerics.iteration` consume them without SN-specific
plumbing.  The within-group inner solve is built once from a single
source of truth — the :func:`~orpheus.sn.coupled_system.build_within_group_system`
builder assembles the :class:`~orpheus.sn.coupled_system.WithinGroupSystem`
record: the posed loss :math:`A` on its carrier space, the bound
LEAVES it is the signed sum of
(:class:`~orpheus.sn.coupled_system.SNLossFactors` — :math:`L+C`,
:math:`S`, :math:`N_{2n}`, :math:`B_a`, plus System B's quartet on a
carrying mesh), and — since 2026-09-13 — the pencil's right-hand side
:attr:`~orpheus.sn.coupled_system.WithinGroupSystem.production`, which
is :math:`F` posed on that same carrier.  :math:`F` is **not** a summand
of :math:`A`: within-group fission is zero, and the production enters as
the :math:`1/k`-scaled outer source.  It rides the record because a
pencil is a *pair* :math:`(A, F)` on one space, and because the factors'
:attr:`~orpheus.sn.coupled_system.SNLossFactors.fission` entry is a
**reference** to the hub's one :math:`F`, never a second mint
(:ref:`sn-one-fission-per-problem`).  Which of the loss leaves is
inverted and
which is lagged is **not** the record's choice: it is a
:class:`~orpheus.sn.splitting.Splitting` — a Strategy VALUE minted from
the record's factors by the one labelling site
:meth:`~orpheus.sn.splitting.Splitting.from_schedule`, whose derived
members :math:`M` (:attr:`~orpheus.sn.splitting.Splitting.implicit`) and
:math:`N` (:attr:`~orpheus.sn.splitting.Splitting.explicit`) the drivers
consume, and whose law :math:`M - N = A` is certifiable per value
against the record it was minted from
(:ref:`sn-splitting-is-a-strategy-value`).  The pieces are handed to the
**variadic** driver
:math:`\text{Driver}(M^{-1},\,*\text{gains})` (Wave O step
O.2a — the transitional :math:`S + B` fold is retired; see
:ref:`bc-extraction-variadic-driver` in :doc:`/theory/foundations/boundary_conditions`).
:func:`_within_group_krylov` wraps the matching
:class:`~orpheus.numerics.iteration.KrylovAcceleration` — and the posed
record is shared verbatim across the eigenvalue source-iteration
inner (:meth:`SNSolver._solve_source_iteration`), the eigenvalue Krylov
inner (:meth:`SNSolver._solve_krylov`), and both fixed-source paths,
each of which mints the splitting VALUE its own schedule calls for.

.. admonition:: Key Facts
   :class: tip

   * The within-group system is built ONCE **per Problem** — it is the
     :func:`~functools.cached_property`
     :attr:`SNProblem.system <orpheus.sn.problem.SNProblem.system>`,
     and that property is the only call site of
     :func:`~orpheus.sn.coupled_system.build_within_group_system`.  It
     carries the posed loss :math:`A`, the bound leaves it is the signed
     sum of (:class:`~orpheus.sn.coupled_system.SNLossFactors`), the
     carrier space, and the posed production
     :attr:`~orpheus.sn.coupled_system.WithinGroupSystem.production`.
     Fission is never inside the swept operator; it enters as the
     :math:`1/k`-scaled outer source.  ``[M]`` the forward eigenvalue
     path built this record once per **outer step** until 2026-09-13 and
     builds it **once** now (:ref:`sn-the-problem-poses-its-pencil`).
   * **The Problem's LAST step is its pencil** (since 2026-09-13 — the
     consumers campaign's step 2, unit C3b).
     :attr:`SNProblem.pencil <orpheus.sn.problem.SNProblem.pencil>`
     is the :class:`~orpheus.numerics.pencil.OperatorPencil`
     :math:`(A, F)` on one space, and
     :attr:`~orpheus.sn.problem.SNProblem.eigen_posing` is the
     :class:`~orpheus.numerics.posing.EigenPosing` over it with the
     :math:`k` map.  The pencil carries **no** inverse and **no**
     resolvent: how it is inverted is the Strategy's
     (:ref:`the-operator-pencil`).  Since 2026-09-14 the solver
     **consumes** the posing:
     :class:`~orpheus.numerics.iteration.KEigenvalue` takes
     ``(posing, implicit, explicit)`` and its three estimators read the
     pencil the posing carries — a principled ULP-level re-baseline of
     the :math:`k` estimator, measured at
     :ref:`sn-the-pencil-reaches-production`.  ⚠ Until 2026-09-14 that
     bullet read *"the solver does not consume the posing yet —
     KEigenvalue still takes the operator triple"*; it was true for one
     commit.
   * **There is ONE** :math:`F` **per Problem, and it lives on the hub**
     (since 2026-09-13 — the consumers campaign's step 2).
     :attr:`SNProblem.fission <orpheus.sn.problem.SNProblem.fission>`
     is the composite :math:`F = \chi\otimes\nu\Sigma_f` on the full
     field; the forward k-outer applies its derived energy face
     ``F.isotropic_energy`` and the adjoint daggers the composite, so
     the two faces cannot disagree.  Before it, the forward and the
     adjoint minted two :math:`F`\ s on two spaces and
     :math:`F_{\rm adjoint}` could not be :math:`F^{\dagger}`
     (:ref:`sn-one-fission-per-problem`).
   * **The splitting** :math:`A = M - N` **is a Strategy VALUE, not a
     member of the posed record** (since 2026-09-13 — the consumers
     campaign's step 2).  Its primitive is the LABELLED TERM SET —
     each leaf carries the sign it has in :math:`A` and exactly one
     label, implicit or explicit — and :math:`M`, :math:`N` are
     *derived* from the labelling, never stored beside it
     (:ref:`sn-splitting-is-a-strategy-value`).  Two labellings ship:
     Jacobi (every geometry) and boundary Gauss-Seidel (multi-D
     Cartesian, seedless), and a value certifies itself against the
     Problem by :meth:`~orpheus.sn.splitting.Splitting.law_residual`.
   * Both inner paths — source iteration and Krylov — consume that SAME
     posed record over the SAME one-walk discretization (matvec ≡
     :term:`sweep`, #206 Phase C): the same solution **set**, different
     rate and memory.  On a closed reflective diamond box that set is a
     *manifold*, not a point — :math:`A` is exactly singular there — so
     the two arms return different **members** and every exit projects
     onto the canonical one (:ref:`sn-loss-kernel-gauge`,
     :ref:`sn-exit-gauge`).
   * :meth:`~orpheus.sn.solver.SNSolver.compute_keff` reports **fission
     production over net removal** (:eq:`sn-keff-update`) — the
     eigenvalue of the map the inner solve actually poses (#291, #259).
     Leakage reads the typed boundary trace; reflective faces are a
     structural zero, so lattice anchors hold to the ULP.
   * Homogenisation and condensation are the frame page's
     Petrov-Galerkin consumers; the SN layer only orchestrates:
     :meth:`Solution.homogenize <orpheus.sn.solution.Solution.homogenize>`
     returns a mesh-coupled ``MaterialMesh``;
     :meth:`Solution.condense <orpheus.sn.solution.Solution.condense>`
     returns a portable ``dict[int, Mixture]`` — the
     condense/homogenize asymmetry law
     (:ref:`sn-condense-homogenize-asymmetry`).


The within-group inner solve consumes the primitives directly
-------------------------------------------------------------

:class:`SNSolver`'s within-group inner solve **is** the
:class:`~orpheus.numerics.iteration.SourceIteration` /
:class:`~orpheus.numerics.iteration.KrylovAcceleration` primitive — not
a verbatim replica of its loop.
:meth:`SNSolver._solve_source_iteration` constructs a
:class:`SourceIteration` from the :func:`~orpheus.sn.coupled_system.build_within_group_system` SSoT and
runs it; :meth:`SNSolver._solve_krylov` constructs a
:class:`KrylovAcceleration` from :func:`_within_group_krylov` and runs
that.  The Layer-3 resolvent of the SN row in the
:ref:`eigenvalue-posing` architecture is exactly these primitive
instances.

The primitive is **type-agnostic and angular-capable**: it operates on
the typed :class:`~orpheus.transport.timed_full_field.TimedFullField`
composite, which carries the full :term:`angular flux` on its bulk.  Pℓ
anisotropic scattering therefore rides the angular bulk with no special
plumbing — :meth:`ScatteringOperator.apply` on the timeless
:class:`~orpheus.transport.full_field.FullField` operator carrier (the
driver's :class:`~orpheus.transport.timed_full_field.TimedFullField` iterate
reaches it via MRO) reads the angular moments off the composite and builds
the anisotropic source inside :meth:`ScatteringOperator.apply
<orpheus.transport.operators.transfer.TransferOperator.apply>` — the
:math:`\ell \ge 1` half being the redistribution body the binding's ends
select (``TransferOperator._redistribute_ordinates`` on the angular end) —
all inside the primitive's normal RHS path.  There is **no scalar-flux
limitation** and **no pending "Approach A" cleanup**: the earlier
framing — that :class:`SourceIteration` carried only :term:`scalar flux` and SN
had to replicate the loop verbatim until the angular state could be
threaded through — was a property of an interim scalar-only carrier
that the typed composite retired.  The
``.claude/skills/algebra-of-record`` "Branch 2 implements the same
operator algebra" discipline is satisfied: SN is the discretized
Branch-2 consumer of the shared primitive, not a parallel loop.

The (L + C − S − N₂ₙ − B)·ψ = (1/k)·F·ψ framing at the solver level
-------------------------------------------------------------------

Beyond driving the within-group inner solve, the :math:`(L+C,\ S,\ F)`
framing organises the solver's outer API surface:

* :meth:`SNSolver.compute_fission_source` returns
  :math:`F\,\phi/k` — a thin delegator to
  ``problem.fission.isotropic_energy.apply`` (the hub's one :math:`F`,
  read at its scalar face) with the :math:`1/k` outer-loop scaling
  applied at the solver level.
* :meth:`SNSolver.solve_fixed_source` solves
  :math:`(L+C-S-N_{2n}-B)\,\psi = q_{\rm ext}`
  (:eq:`sn-within-group-with-n2n`; with :math:`q_{\rm ext}` the
  fission source built by ``compute_fission_source``).  Two paths:

  * ``inner_solver="source_iteration"`` — sweep-driven fixed-point
    iteration; the resolvent :math:`(L+C)^{-1}` is realised by the
    one-walk WDD sweep.
  * ``inner_solver="krylov"`` — GMRES on the honest within-group
    matvec, with the sweep resolvent as preconditioner — the same
    one-walk discretization either way (matvec ≡ sweep, #206 Phase C).

* :meth:`SNSolver.compute_keff` returns **fission production over net
  removal**, :eq:`sn-keff-update` — the volume-weighted method-layer
  functional :math:`R_{\nu\Sigma_f}(\phi) / (R_{\Sigma_a}(\phi) + L -
  E_{2n}(\phi))`, derived in :ref:`sn-keff-estimator` below.  The
  SN-specific volume weighting lives here (in the typed
  :class:`~orpheus.transport.reaction_rate_functional.IntegratedReactionRate`
  fields) — one honest realization of the same discipline the
  operator-form :meth:`KEigenvalue.compute_keff` spells with the
  measure absorbed into the operators' action.  (Pre-#291 this method
  returned the leakage-blind :math:`\sum F\phi V / \sum \Sigma_a\phi V`
  ratio; see :ref:`sn-keff-estimator` for why that was a
  non-eigenvalue on any vacuum-bounded problem.)

The solver-level :math:`1/k` scaling (in
:meth:`~SNSolver.compute_fission_source`) and the volume-weighted
eigenvalue estimate (in :meth:`~SNSolver.compute_keff`) are exactly the
points where SN's specifics live; the rest of the solver is
operator-algebra coordination over the canonical
:func:`~orpheus.numerics.eigenvalue.power_iteration` boundary.  These
two K-specific hooks are also precisely *why* the Layer-4 loop is not
yet literally K/α-agnostic — relocating the eigenvalue scaling into the
algorithm is the first step of the α-wave (see the honest-scope caveat
in :ref:`eigenvalue-posing`).

The eigenvalue :math:`\keff` is determined by **power iteration**: an
outer loop updates :math:`k` from the net-removal balance
:eq:`sn-keff-update` (fission production over absorption + leakage −
:math:`(n,2n)` emission), with an inner loop that solves the
within-group scattering problem.

.. _sn-keff-estimator:

The reported eigenvalue: fission production over net removal
------------------------------------------------------------

:meth:`~orpheus.sn.solver.SNSolver.compute_keff` reports the eigenvalue
of the problem the inner solve **actually poses**.  This is the SN
symptom (#291) and the MoC/CP/homogeneous root (#259) of a single
discipline: *the reported* :math:`k` *must be the eigenvalue of the
fixed-source map every method scales only fission by* :math:`1/k`
*through* — scattering and the :math:`(n,2n)` emission are plain gains
assembled **inside** :meth:`~orpheus.sn.solver.SNSolver.solve_fixed_source`,
never scaled by :math:`1/k`.  An estimator that disagrees with its own
posed problem converges cleanly and silently to a **non-eigenvalue
ratio**.

.. math::
   :label: sn-keff-update

   k \;=\; \frac{R_{\nu\Sigma_f}(\phi)}
                {R_{\Sigma_a}(\phi) \;+\; L \;-\; E_{2n}(\phi)}

.. (V&V scope note) Governing/definitional identity: the reported k
   IS the eigenvalue of the posed fixed-source map, not a solver
   eigenvalue-correctness claim against an external analytical reference
   (that rests on the multi-group heterogeneous L1/L2 references
   elsewhere on this page). The label is wired to the cross-engine
   consistency gate tests/sn/eigenvalue/test_keff_estimator_gate.py
   (reported k == the converged fixed-point map ratio k* = P(Mφ*)/P(φ*),
   map-ratio ground-truth noise ≤ 2e-11) with in-file mutation teeth.

The three terms are typed volume-integrated reaction-rate functionals
and one boundary functional:

* **Numerator** :math:`R_{\nu\Sigma_f}(\phi) = \int_V \nu\Sigma_f\,\phi\,dV`
  — the fission production, the typed
  :class:`~orpheus.transport.reaction_rate_functional.IntegratedReactionRate`
  over :math:`\nu\Sigma_f` (the :math:`\phi^\dagger\!=\!1` degenerate of
  the homogenization Petrov–Galerkin bilinear).  The :math:`(n,2n)`
  emission is **not** production here — contrast
  :meth:`~orpheus.sn.solver.SNSolver.compute_production_rate`, the
  ERR-052 renormalisation scale anchor, which keeps *total* physical
  production (fission **plus** the :math:`(n,2n)` emission).  The role
  split is the load-bearing #259 correction: the same physical
  :math:`(n,2n)` neutrons are a **scale** quantity for the outer
  renormalisation but a **removal-reduction** in the eigenvalue balance.
* **Absorption** :math:`R_{\Sigma_a}(\phi) = \int_V \Sigma_a\,\phi\,dV`,
  with :math:`\Sigma_a = \Sigma_f + \Sigma_c + \Sigma_L +
  \sum_{g'}\Sigma_{2,g\to g'}` — i.e. ``absorption_xs`` counts the
  :math:`(n,2n)` **collision once** (the neutron is removed from its
  incident group by the collision).  See
  :attr:`~orpheus.data.macro_xs.mixture.Mixture.absorption_xs`.
* **Leakage** :math:`L` — the net vacuum-boundary outflow (below).  On a
  reflective (lattice) problem it is a **structural zero**.
* **Emission** :math:`E_{2n}(\phi) = \int_V \sum_{g,g'} 2\,\Sigma_{2,g'\to
  g}\,\phi_{g'}\,dV` — the :math:`(n,2n)` **emission** (two neutrons out
  per collision; the factor 2).  A gain, so it **reduces** net removal.

The net :math:`(n,2n)` effect on removal is therefore
:math:`\underbrace{\sum_{g'}\Sigma_{2,g\to g'}}_{\text{collision, in }\Sigma_a}
- \underbrace{2\Sigma_2}_{E_{2n}} = -\Sigma_2` — **one extra neutron
gained** per collision, exactly the physics of a neutron-doubling
reaction.

**Balance identity (divergence-telescoping).**  The angle- and
group-summed discrete cell balance for cell :math:`i` in the posed
eigenproblem is

.. math::
   :label: sn-keff-cell-balance

   \underbrace{\sum_{f\in\partial i}\!\bigl(\textstyle\sum_g J_g\bigr)\,\Delta A_f}
              _{\text{net face flow}}
   \;+\; \Sigma_{t,i}\,\phi_i\,V_i
   \;=\; \frac{1}{k}\,R_{\nu\Sigma_f,i}
        \;+\; \Sigma_{s,i}\,\phi_i\,V_i
        \;+\; E_{2n,i}

.. (vv-status rationale) Derivation step (the divergence-telescoping cell
   balance). Its terminal result sn-keff-update is verified by the k* map-ratio
   gate (tests/sn/eigenvalue/test_keff_estimator_gate.py); definitional.
.. vv-status: sn-keff-cell-balance documented

(streaming + total collision on the left; the isotropic fission source
scaled by :math:`1/k`, plus the *unscaled* scatter and :math:`(n,2n)`
gains, on the right).  Summing over **all** cells, every interior face
is shared by two cells with opposite outward normals and equal current
(continuity), so the interior face-flow terms **telescope to zero** —
only the domain-boundary faces survive, and their sum is the net
leakage :math:`L`.  With :math:`\Sigma_t - \Sigma_s = \Sigma_a` this
collapses to

.. math::

   \frac{R_{\nu\Sigma_f}(\phi)}{k} \;=\; R_{\Sigma_a}(\phi) \;+\; L
                                        \;-\; E_{2n}(\phi),

which is :eq:`sn-keff-update` rearranged.  This is the same discrete
divergence discipline the diffusion page states as
:math:`\mathbf 1^{\mathsf T}(C-S)=\Sigma_a` with interior-face
telescoping (see :ref:`diffusion-leakage-boundary-leaves`); SN and
diffusion report the *same* balance-law eigenvalue, differing only in
how the streaming operator is discretised.

The leakage functional
~~~~~~~~~~~~~~~~~~~~~~~~

.. math::
   :label: sn-leakage-functional

   L \;=\; \sum_{f\,\in\,\text{vacuum}} \oint_{f} dA\,
           \sum_g J_g(\mathbf{r}_f)\,,
   \qquad
   J_g \;=\; \sum_m (\Omega_m\cdot\hat n_f)\, w_m\, \psi_{m,g}

is the face-area integral of the boundary trace's **net outward
current**.  The angular-to-scalar reduction :math:`J_g` is
:meth:`AngularBoundaryFlux.net_current
<orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux.net_current>`
— the single source of the :math:`\Omega\cdot\hat n\,w` contraction, the
angular sibling of the scalar trace's :math:`J = J^+ - J^-`
(:meth:`ScalarBoundaryFlux.net_current
<orpheus.transport.fields.scalar_boundary_flux.ScalarBoundaryFlux.net_current>`).
It is spelled through the trace space's own atoms — the signed
projection table
:attr:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.omega_dot_n`
and the :math:`|\Omega\cdot\hat n|\odot w` partial-current metric (using
the identity :math:`\operatorname{sign}(\Omega\cdot\hat n)\cdot
|\Omega\cdot\hat n|\,w = \Omega\cdot\hat n\,w`) — so no consumer
re-derives the cosine weighting.  Tangential :term:`ordinates <ordinate>` carry zero weight
and drop out.

The face measure :math:`dA` is supplied by
:meth:`SNSolver._face_area_of`, matching the cell
``volume_measure`` exactly so the balance identity closes:

.. list-table:: Boundary-face measure by geometry
   :header-rows: 1
   :widths: 30 30 40

   * - Geometry
     - Face measure :math:`\Delta A`
     - Source
   * - 1-D slab
     - :math:`1` (per unit cross-section)
     - :attr:`MaterialMesh.areas <orpheus.transport.mesh.material_mesh.MaterialMesh.areas>`
   * - 1-D cylinder
     - :math:`2\pi R` (per unit height)
     - ``MaterialMesh.areas``
   * - 1-D sphere
     - :math:`4\pi R^2`
     - ``MaterialMesh.areas``
   * - 2-D Cartesian
     - transverse edge-cell widths (unit depth)
     - ``mesh.axes`` transverse extent
   * - 3-D Cartesian
     - :math:`\Delta A_{\mathbf c} = \prod_{j\ne a}\Delta_j[c_j]`
       (transverse-area outer product)
     - ``mesh.axes`` transverse extents

The :math:`d \ge 2` Cartesian arms are ONE generic body: the outer
product of the *other* axes' edge widths in **ascending axis order** —
the same codimension-1 enumeration as
:func:`~orpheus.transport.mesh.axis.face_shape`, so the measure array
broadcasts cell-for-cell against the ``(ng, *face_spatial)`` net
current, and the 2-D width vector is just the single-transverse-axis
degenerate (bit-identical to the pre-3-D spelling).

The 3-D arm originally shipped as a **typed refusal**
(``NotImplementedError``): guessing the transverse product's cell
ordering would silently mis-weight the leakage sum, and Cardinal Rule 1
forbids returning a wrong-but-plausible number.  The wire landed
2026-07-13 when the first 3-D vacuum eigenvalue consumer arrived (the
d=3 Mode-9 G-S≡Jacobi gate), with the ordering pinned twice in
``tests/sn/eigenvalue/test_keff_estimator_gate.py``: an **object-level
pin** (face measure ≡ the boundary layer's ``volumes / Δ_axis``, the
mesh's own ascending-axis enumeration — vv Mode-12 discipline: pin the
object, not only the k functional) and the **k* map-ratio gate** on a
Mode-2 asymmetric all-vacuum box, whose teeth are proven by permanent
in-process mutants — a reversed transverse enumeration moves the
reported k by a measured **13.9 %** against the estimator-independent
:math:`k^*` (clean agreement :math:`6\times10^{-10}`), and a transposed
enumeration crash-REDs on the broadcast.  A trace carrying a
``#251`` transverse face-moment tail is refused loudly at the
consumption site (the face integral must consume ONLY the
transverse-average moment — higher Legendre face moments integrate to
zero over each face cell — and that slot-0 read has no consumer yet).

Reflective faces are a structural zero
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The reflective law equates a face's inflow to its reflected outflow
**exactly**, so the net current there vanishes *by construction*.
:meth:`~orpheus.sn.solver.SNSolver._boundary_leakage_rate` therefore
**skips** reflective faces — it never accumulates them, rather than
accumulating a value that ought to be zero but carries
:math:`\pm`-cancelling angular-sum floating-point noise.

This is a deliberate design choice with a bit-level payoff.  On an
all-reflective (lattice) problem :math:`L` is a structural ``0.0``, and
on a :math:`\Sigma_2`-free mixture :math:`E_{2n}` is exactly ``0.0`` (the
per-material :math:`(n,2n)` loop adds nothing), so
:eq:`sn-keff-update` reduces **bit-identically** to the historical
lattice functional ``production / absorption``.  Every pre-existing
reflective eigenvalue anchor is preserved to the last ULP — the
unification adds terms that vanish structurally, not numerically, on the
cases it must not perturb.

The scale bridge: trace of the last inner solve
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The leakage term reads the **typed** boundary trace of the last inner
solve (``self._inner.iterate.boundary`` — the
:class:`~orpheus.sn.solver.InnerSolve` record each inner leaves behind),
whereas the numerator/denominator
reaction rates consume the bare-array flux :math:`\phi` the estimator is
handed.  These two representations can be at **different scales**:
:func:`~orpheus.numerics.eigenvalue.power_iteration` renormalises
:math:`\phi` to unit production rate *between* the inner solve and the
:math:`k`-update (ERR-052), so the stored trace belongs to the
**un-renormalised** last iterate while the estimator's :math:`\phi` is
its renormalised multiple.

Leakage is degree-1 homogeneous in :math:`\psi`, so the fix is a single
rescale by the fission-production ratio of the two fluxes
(``self._phi_of_trace``, stored alongside the trace at **both**
inner-path returns) — exactly ``1.0`` when the caller passes the
returned flux itself.  The **contract** is therefore: the flux handed to
:meth:`~orpheus.sn.solver.SNSolver.compute_keff` must be a scalar
multiple of the last inner solve's flux (true for ``power_iteration`` and
for every manual solve-then-estimate loop).

If a vacuum face exists but **no** inner solve has stored a trace,
:meth:`~orpheus.sn.solver.SNSolver._boundary_leakage_rate` raises
``RuntimeError`` — the leakage cannot be answered honestly, and silently
returning it as zero would *reproduce the #291 omission*.  Fail loud;
never return a non-eigenvalue.

The R7 :math:`(n,2n)` convention fork
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The historical spelling put the :math:`(n,2n)` emission in the
**numerator** as production,

.. math::
   :label: sn-keff-old-n2n

   k_{\text{old}} \;=\; \frac{R_{\nu\Sigma_f} + E_{2n}}{R_{\Sigma_a}},

.. (vv-status rationale) Definitional/historical contrast: the superseded
   (n,2n)-in-numerator spelling. No current code implements it; its bias is
   characterised in the pre-fix table (#291 commit d1daaac).
.. vv-status: sn-keff-old-n2n documented

which is a **non-eigenvalue** of the posed map whenever
:math:`\Sigma_2 \neq 0` *and* :math:`k \neq 1`.  The reason is exactly
the posing asymmetry: the inner solve scales **only** fission by
:math:`1/k`; the :math:`(n,2n)` emission is an *unscaled* gain in the
sweep source.  So the eigenvalue of that map is
:math:`k^\star = R_{\nu\Sigma_f}/(R_{\Sigma_a} - E_{2n})` (reflective,
:math:`L=0`), and putting the *unscaled* emission in the numerator does
not recover it.  Writing :math:`f = R_{\nu\Sigma_f}`,
:math:`a = R_{\Sigma_a}`, :math:`e = E_{2n} = 2s_2` and substituting
:math:`f = k^\star (a - e)`:

.. math::
   :label: sn-keff-old-bias

   k_{\text{old}}
   \;=\; \frac{k^\star (a - e) + e}{a}
   \;=\; k^\star \;+\; \frac{2 s_2\,(1 - k^\star)}{a}.

.. (vv-status rationale) Mathematical identity (the derived bias of the
   superseded estimator, k_old - k* = 2 s_2 (1 - k*) / a). Historical
   characterisation; no current implementing code.
.. vv-status: sn-keff-old-bias documented

The two agree only when :math:`s_2 = 0` (no :math:`(n,2n)`) or
:math:`k^\star = 1` (critical).  For a supercritical
:math:`k^\star > 1` the correction is negative
(:math:`k_{\text{old}} < k^\star`); for subcritical, positive.  The MoC
and CP pages carry the same fork
(:eq:`moc-keff-update`, :eq:`cp-keff-update`); CP was the one member
already spelled on net removal.

What was tried and found
~~~~~~~~~~~~~~~~~~~~~~~~~~

The #291 bias was characterised pre-fix (commit ``d1daaac``, Gauss–
Legendre :math:`n=8`, map-ratio ground truth noise :math:`\le 2\times
10^{-11}`) across the five gate configurations:

.. list-table:: Pre-fix reported :math:`k` vs the posed-problem eigenvalue :math:`k^\star`
   :header-rows: 1
   :widths: 40 18 18 24

   * - Configuration
     - Pre-fix reported
     - Posed :math:`k^\star`
     - Bias
   * - homog. 2G vacuum slab (width 8)
     - 1.83767525
     - 0.98163269
     - :math:`+87.2\%` (:math:`L/A = 0.872`)
   * - het. vacuum sphere P\ :sub:`0`
     - 0.86484694
     - 0.70601977
     - :math:`+22.5\%`
   * - het. vacuum sphere P\ :sub:`1`
     - 0.85080423
     - 0.67876772
     - :math:`+25.3\%`
   * - reflective control (:math:`\Sigma_2=0`)
     - 1.87500000
     - 1.87500000
     - :math:`\equiv` (bias :math:`1.2\times 10^{-10}`)
   * - reflective :math:`\Sigma_2\neq 0`
     - 1.92857143
     - 2.61278195
     - :math:`-26.2\%` (R7 defect)

The two failure classes are visible in one table: the vacuum rows are
the **leakage omission** (#291) — the reported :math:`k` overshoots by
the leakage-to-absorption ratio :math:`L/A`; the last row is the **R7
:math:`(n,2n)` convention** — zero leakage, yet a
:math:`-26.2\%` error because the emission was mis-posed as production.
The reflective-control row is exactly the bit-identity guarantee above.
The exact check on the R7 row is
:math:`0.78/(0.5185 - 0.2200) = 2.61278`, and
:math:`(0.78 + 0.2200)/0.5185 = 1.92857` reproduces the old value —
matching :math:`k_{\text{old}} = k^\star + 2s_2(1-k^\star)/a` term for
term.

Post-fix, reported :math:`k` and the map-ratio :math:`k^\star` agree to
:math:`\le 6\times 10^{-10}` on all five.  The P\ :sub:`0`\ –P\ :sub:`1`
sphere gap :math:`\Delta` roughly **doubled** (:math:`1.404\times
10^{-2} \to 2.725\times 10^{-2}`) but stays inside the diagnostic
:math:`(10^{-3}, 5\times 10^{-2})` band — the P\ :sub:`1` anisotropic
correction is now measured against the correct eigenvalue on both
solves.

The V&V decision was a **principled re-baseline** (per ``vv-principles``
bit-identity-vs-principled-equivalence): the old reported :math:`k` was a
*different functional* from the posed problem's eigenvalue, so the new
value is not a regression to be tolerance-matched but a correction to be
verified against a structurally-independent reference (the fixed-point
map ratio).

Verification
~~~~~~~~~~~~

The permanent gate is
``tests/sn/eigenvalue/test_keff_estimator_gate.py``: it asserts the
reported :math:`k` equals the converged fixed-point map ratio
:math:`k^\star = P(M\phi^\star)/P(\phi^\star)` across the four physics
regimes — {vacuum slab, vacuum sphere (pinning the :math:`4\pi R^2`
face-area convention), reflective bitwise-degenerate, reflective
:math:`\Sigma_2\neq 0`} — with **in-file mutation teeth**: a
leakage-drop mutation reds the vacuum legs while staying bitwise-green
on reflective; a leakage sign-flip crash-reds through the scale-bridge
guard; and the old :math:`(n,2n)`-in-numerator convention reds the
:math:`\Sigma_2\neq 0` leg.

This is a **consistency** gate: the map ratio is the structurally-
independent ground truth for "does the estimator return the eigenvalue
of its own posed map", and is blind by construction to *which*
eigenvalue that is.  The SN solver's eigenvalue **correctness** — that
the posed map's eigenvalue is the *physically right* :math:`k` — rests
on the multi-group heterogeneous L1/L2 references in
:doc:`/theory/verification/sn`, not on this gate.

Two Inner Solvers
-----------------

**Source iteration (sweep-based):**

- Operator: :math:`(L+C)^{-1}` (the one-walk WDD sweep)
- Iterate: the typed field composite (angular bulk + boundary trace)
- Fixed-point: :math:`\psi^{(k+1)} = (L+C)^{-1}(S\,\psi^{(k)} +
  N_{2n}\,\psi^{(k)} + B\,\psi^{(k)} + q_{\rm ext})` — the Jacobi
  value's :attr:`~orpheus.sn.splitting.Splitting.explicit` triple
  ``(S, N2N, B_a)``
- Convergence rate: spectral radius of
  :math:`(L+C)^{-1}(S+N_{2n}+B)` — the
  :term:`scattering ratio` :math:`c` (:doc:`slab_one_group`)
- Cost per iteration: one transport sweep
- Works for all geometries

**Krylov (direct operator):**

- Operator: the honest :math:`(L+C-S-N_{2n}-B)` applied matrix-free (its
  :math:`(L+C)` piece via :meth:`StreamingCollisionOperator.apply` — the same
  one-walk discretization the sweep realises; L21 matvec ≡ sweep)
- Iterate: the same typed composite; GMRES additionally stores its
  Krylov basis (``restart`` × the composite's ``n_dof`` — the ERR-053
  sizing family)
- System: :math:`(L+C-S-N_{2n}-B)\,\psi = q_{\rm ext}` — scattering,
  the :math:`(n,2n)` gain and the
  boundary gain live in the operator, not the lagged source;
  :math:`q_{\rm ext}` is the :math:`1/k`-scaled fission source
- Convergence: GMRES with sweep preconditioner, typically ~100
  Krylov iterations at ``tol=1e-4`` (always converges)
- Available for all geometries (Cartesian, spherical, cylindrical)

Wave E Round 2 (Issue #164) replaced the legacy BiCGSTAB FD-operator
path with this Krylov path.  See the Krylov alternative in
:doc:`slab_one_group` for the full discussion.

The two paths share the **one** loss-representation discretization
(matvec ≡ sweep, #206 Phase C), so they solve the same system and
converge to the same **solution set**; the Wave-D-era design in which
they carried different spatial closures — and disagreed on coarse-mesh
:math:`\keff` — is recorded in the two-closure history
(:ref:`loss-rep-history`).

.. note::

   ⛔ **That read "the same fixed point" until 2026-08-15 (#344), and a
   set is not a point.**  On a closed reflective diamond box
   :math:`A = L+C-S-N_{2n}-B` is **exactly singular**, so the two arms
   legitimately return different members of a solution manifold:
   ``[M]`` on an all-reflective 2-D absorber box with a uniform
   isotropic source, source iteration under boundary Gauss-Seidel
   returns a trace carrying :math:`6.08\times10^{-2}` of
   :math:`\ker A` while Krylov carries :math:`4.1\times10^{-14}`, and
   the difference between them lies **entirely** in :math:`\ker A`
   (:math:`\lVert\Pi d\rVert/\lVert d\rVert = 1.000000`).  The **bulk**
   really is arm-invariant — the kernel is pure-trace — and every entry
   now projects the returned trace onto the canonical member, so the
   sentence is true again of what a caller receives.  Derivation:
   :ref:`sn-loss-kernel-gauge`; exit behaviour: :ref:`sn-exit-gauge`.

.. _sn-splitting-is-a-strategy-value:

The splitting ``A = M - N`` is a Strategy VALUE
-----------------------------------------------

The two inner solvers above do not consume the posed loss :math:`A`
raw.  They consume a **splitting** of it — a decomposition
:math:`A = M - N` in which :math:`M` is inverted every step and
:math:`N` is evaluated on the previous iterate (Hackbusch 2016, §11).
Source iteration *is* that decomposition
(:math:`\psi \leftarrow M^{-1}(q + \sum_i N_i\psi)`); Krylov uses
:math:`M^{-1}` as its preconditioner and applies :math:`M - N`
matrix-free.

Choosing which leaf is inverted and which is lagged is a **solver**
decision, not a property of the equation.  Since 2026-09-13 (the
consumers campaign's step 2) the codebase says so structurally: the
Problem's posed record
:class:`~orpheus.sn.coupled_system.WithinGroupSystem` carries the loss,
its carrier space, and the bound leaves
(:class:`~orpheus.sn.coupled_system.SNLossFactors`) — and *nothing
else*.  The splitting is a separate frozen value,
:class:`~orpheus.sn.splitting.Splitting`, minted from those leaves by
one labelling site,
:meth:`~orpheus.sn.splitting.Splitting.from_schedule`.

.. admonition:: Key facts for this section
   :class: tip

   * The **primitive is the labelling**, not the products.  A
     :class:`~orpheus.sn.splitting.LossTerm` is an operator together
     with the coefficient :math:`\pm 1` it carries in :math:`A`; a
     splitting assigns each term exactly one label, implicit or
     explicit, and :math:`M` and :math:`N` are *derived*
     (:eq:`sn-splitting-labelled-terms`).  Nothing stores :math:`M`
     and :math:`N` beside the labelling, so they cannot disagree with
     it.
   * **Two labellings ship**: Jacobi (every geometry, and the only one
     admitted on a seed-carrying mesh) and boundary Gauss-Seidel
     (multi-D Cartesian, seedless).  They are two splittings of the
     *same* :math:`A`, so they share a fixed-point *set* and differ
     only in rate — the pair is what a Mode-9 invariance gate needs.
   * The law :math:`M - N = A` is a property of the value and is
     checkable per value:
     :meth:`~orpheus.sn.splitting.Splitting.law_residual`.  It is
     **bit-exact** on the seedless arm (the Gauss-Seidel split writes
     disjoint rows, so no addition is reordered) and round-off-close on
     the carrying arm (the block grid re-associates).
   * ⚠ It is a splitting, **not a regular splitting** in Varga's sense,
     so no comparison theorem bounds its rate and boundary
     Gauss-Seidel is measurably slower than Jacobi on some meshes:
     :ref:`sn-boundary-gs-not-regular`.

The defect this retired: a posing that advertised a splitting it did not use
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

From 2026-07-28 until 2026-09-13 the posed record carried two further
members, ``implicit_operator`` (:math:`M`) and ``explicit_gains``
(:math:`N`), assembled by the builder.  That was wrong twice over.

It was wrong **in principle** because a posing has no schedule: the same
:math:`A` is posed once and then solved by whichever strategy the caller
selected, so a record that names one :math:`M` is answering a question
it was never asked.

It was wrong **in fact**, and that is the part a gate could see.  The
seedless source-iteration driver did *not* run the record's pair.  On a
multi-D Cartesian mesh with ``inner_schedule="gauss_seidel"`` it called a
solver-private helper, ``_select_si_splitting``, which re-derived a
*second* splitting behind the record's back: it split :math:`B_a` under
the octant-group schedule and folded the strictly-lower half into the
implicit operator.  So the record advertised
:math:`M = (L+C)`, :math:`N = (S, N_{2n}, B_a)` while the driver ran
:math:`M = (L+C) - B_{\rm lower}`, :math:`N = (S, N_{2n}, B_{\rm upper})`
— a twin, in exactly the Cardinal-Rule-2 sense, and the record's claim
was simply false for one of the two shipped schedules.

That twin was tracked as **R7 of the operator/strategy campaign** (not
to be confused with the ``(n,2n)`` R7 of :ref:`sn-keff-estimator`, nor
with the #310 R7 schedule-reverse transpose).  It was pinned by a
``xfail(strict=True)`` row in
``tests/sn/architecture/test_stage_separation.py`` so the fix could not
land silently; step 2 made the row XPASS, and the marker was deleted the
same day.  What keeps teeth now is the law itself, below, plus the
positive rows that assert the driver consumes the value's own operators
by identity.

.. note::

   The record's ``implicit_operator`` field had an earlier name,
   ``resolvent``, retired 2026-07-28 as a misnomer — it held the
   *un-inverted forward* :math:`M`, whereas a resolvent is
   inverse-like.  The two renames are one story: first the field was
   named honestly, then it was found not to belong on that object at
   all.  See the crosswalk row in
   :doc:`/theory/conventions/notation`.

The labelled term is the primitive
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The within-group loss on the loss-sign convention is the signed sum
:eq:`sn-within-group-with-n2n`, plus System B's quartet on a
seed-carrying mesh (:ref:`coupled-block-operator`):

.. math::

   A \;=\; \underbrace{(L+C)}_{+1} \;-\; S \;-\; N_{2n} \;-\; B_a
   \qquad \bigl(
   +\,A_{AB}\ \text{(seeding)},\;
   -\,E\ \text{(emission)},\;
   +\,A_{BB}\ \text{(march)},\;
   -\,B_b
   \bigr).

A **term** is therefore an operator together with the coefficient it
carries here, and that pair is the object
:class:`~orpheus.sn.splitting.LossTerm` reifies: ``LossTerm(operator,
sign)`` with ``sign`` in :math:`\{+1, -1\}` (the constructor refuses
anything else, so a corrupt coefficient is unspellable).  A splitting is
then a **partition of the term set into two labels**, implicit
:math:`\mathcal I` and explicit :math:`\mathcal E`, and its two members
are derived from that partition:

.. math::
   :label: sn-splitting-labelled-terms

   A \;=\; \sum_{k} s_k\,T_k ,
   \qquad
   M \;=\; \sum_{k \in \mathcal I} s_k\,T_k ,
   \qquad
   N \;=\; -\!\!\sum_{k \in \mathcal E} s_k\,T_k ,
   \qquad
   \mathcal I \sqcup \mathcal E = \{k\} .

.. (vv-status rationale) Representational identity: it states how the
   shipped ``Splitting`` value derives its two members from the labelling
   it holds, which is a definition rather than a solver claim.  The
   verifiable content is the law :eq:`sn-splitting-law` it makes
   immediate, gated by
   ``tests/sn/solve/test_boundary_gs_is_a_coherent_splitting.py::test_both_schedules_are_splittings_of_the_SAME_A``
   (exhaustive dense assembly, both schedules) and
   ``tests/sn/architecture/test_stage_separation.py::test_reconstruction_identity_A_equals_M_minus_N``
   (both arms), with their mutation rows.
.. vv-status: sn-splitting-labelled-terms documented

The minus sign on :math:`N` is what makes the lagged terms **gains**.  A
gain carries :math:`s_k = -1` in :math:`A`; flipping it puts the term on
the right-hand side with a plus, which is what the driver's update
:math:`\psi \leftarrow M^{-1}(q + \sum_i N_i \psi)` needs.  The law then
holds by construction rather than by convention — write it out:

.. math::
   :label: sn-splitting-law

   M - N
   \;=\; \sum_{k \in \mathcal I} s_k T_k
       \;+\; \sum_{k \in \mathcal E} s_k T_k
   \;=\; \sum_{k} s_k T_k
   \;=\; A .

.. (vv-status rationale) Structural identity: ``A = M - N`` is immediate
   from :eq:`sn-splitting-labelled-terms` once the labels partition the
   term set, so the equation states a property of the shipped
   construction, not a claim about a solve.  Its verifiable content —
   that the objects the value actually derives satisfy it on real
   operands — is what
   :meth:`~orpheus.sn.splitting.Splitting.law_residual` computes and
   what the two gate modules named above assert, bit-exactly on the
   seedless arm and at ``nulp=8`` on the carrying one.
.. vv-status: sn-splitting-law documented

Three design consequences follow from carrying the sign on the **term**
rather than on the label, and each is the reason a plausible-looking
alternative was rejected.

*The sign survives a label move.*  Moving a term from
:math:`\mathcal I` to :math:`\mathcal E` changes which sum it appears
in, not its coefficient — so a *transfer* is law-invariant by
construction.  ``[M]`` transferring :math:`B_{\rm lower}` from implicit
to explicit leaves the law residual at exactly ``0.0``.  Had the sign
been positional (implicit terms added, explicit terms subtracted), the
same transfer would have moved the residual by :math:`2\lvert T_k x
\rvert` — i.e. the operation that a later ρ-policy phase needs as its
elementary move would have been unsafe by default.

*The derivation stays on the domain's own algebra.*  A term is stored as
``(operator, sign)`` rather than pre-wrapped in a
:class:`~orpheus.numerics.operator.ScaledOperator` with coefficient
:math:`-1`, because the fold must be able to spell
:math:`(L+C) - B_{\rm lower}` as a *subtraction*.  That expression
dispatches, through
:class:`~orpheus.sn.operators.streaming.StreamingCollisionOperator`'s
own ``__sub__``, to the sweep-invertible scheduled composite whose
``solve`` is the octant-group forward substitution.  An equivalent
``+ (-B_lower)`` would compose a generic
:class:`~orpheus.numerics.operator.OperatorSum` instead, and the result
would not be invertible by a sweep.

*Shape symmetry holds where it must.*  The design rule "a block cannot
move from :math:`M` to :math:`N` if :math:`M` and :math:`N` have
different shapes" is false of the *products* — :math:`M` is one operator
and :math:`N` is a tuple of pieces on the seedless arm — and true of the
*terms*, which is the level the rule was always about.  Because the
primitive is the term set, a transfer has an operand.

The two shipped labellings
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 16 30 30 24

   * - labelling
     - implicit :math:`\mathcal I`
     - explicit :math:`\mathcal E`
     - admitted on
   * - **Jacobi**
     - :math:`L+C` (and, carrying, :math:`A_{AB}` seeding and
       :math:`A_{BB}` march)
     - :math:`S`, :math:`N_{2n}`, :math:`B_a` (and, carrying,
       :math:`E` emission and :math:`B_b`)
     - every geometry; the **only** labelling on a seed-carrying mesh
   * - **Boundary Gauss-Seidel**
     - :math:`L+C`, :math:`B_{\rm lower}`
     - :math:`S`, :math:`N_{2n}`, :math:`B_{\rm upper}`
     - multi-D Cartesian, seedless

Under boundary Gauss-Seidel the *only* operator that splits is
:math:`B_a`.  The octant-group schedule induces a partition of its rows
into a strictly-lower half (inflow rows whose reflected source is
already available within the sweep) and an upper half, via
:meth:`~orpheus.sn.operators.boundary.SNBoundaryOperator.split`, which
returns the named pair :math:`B_{\rm lower} + B_{\rm upper} = B_a` on
**disjoint** rows.  The lower half then moves across the :math:`M/N`
boundary and is absorbed into the implicit part; the collision gains
:math:`S` and :math:`N_{2n}` are lagged under *both* labellings,
because the sweep never re-scatters mid-sweep.

Two refusals make the restriction structural rather than conventional:

* The schedule string becomes a schedule object at exactly one site,
  :func:`~orpheus.sn.splitting.resolve_schedule`, and that site carries
  the geometry gate.  ``"gauss_seidel"`` yields a sequenced schedule
  only when the Problem ``is_cartesian and not is_1d``; on 1-D or
  curvilinear geometries it falls back to Jacobi.  The gate reads the
  genuine condition: before #225 C5.4 it read the proxy ``reduced is
  None``, which happened to be equivalent on 2-D Cartesian and was not
  the property being tested.
* Handing a sequenced schedule to a carrying system
  :meth:`~orpheus.sn.splitting.Splitting.from_schedule` **raises** rather
  than silently re-labelling.  A carrying Problem's boundary is the
  :math:`B_a + B_b` composite, and the ruling that the octant grading
  lives on :math:`B_a` and only there would be violated by splitting it.
  ``resolve_schedule`` never produces that combination, so the refusal
  exists for a direct caller that bypasses it.

Where :math:`M` and :math:`N` come from
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The derivation has two arms, and they differ in *mechanism*, not in the
definition above.

**Seedless — a signed fold.**  The implicit terms are folded left to
right through the operator algebra: the first term *starts* the sum, and
every further term is added or subtracted according to its sign.  Three
properties of that choice matter:

* A single-term implicit set folds to **that operator itself**, not to a
  one-element sum wrapping it.  Identity is preserved, which is what lets
  the stage-separation gates assert that the driver inverted the record's
  own :math:`L+C` object by ``is``, rather than by value.
* A two-term implicit set :math:`\{L+C,\; -B_{\rm lower}\}` folds
  through the subtraction dispatch described above and lands on the
  scheduled composite, so :math:`M^{-1}` is the octant-group forward
  substitution with no branch anywhere in the driver.
* The fold refuses an empty implicit set, and refuses to *open* with a
  gain — the first term must be a :math:`+` term, because the sweepable
  structure is what makes :math:`M` invertible at all.

**Carrying — placement by ends.**  On a seed-carrying mesh every term is
placed into the :math:`2\times2` System-A :math:`\oplus` System-B block
grid **by its own ends**: a term's codomain names its block row and its
domain names its block column, resolved against the members of the
coupled carrier space.  No term carries a block index; the coupled
space's members *are* the addresses.  Terms landing in the same slot are
summed in labelling order, and a term whose end matches zero members —
or more than one — raises, so an unbound end is caught as the posing
defect it is rather than becoming a placement choice.

That yields the honest upper-triangular implicit grid and the gain grid

.. math::

   M \;=\;
   \begin{bmatrix} L+C & A_{AB} \\ \varnothing & A_{BB} \end{bmatrix} ,
   \qquad
   N \;=\;
   \begin{bmatrix} S + N_{2n} + B_a & \varnothing \\ E & B_b
   \end{bmatrix} ,

in which the :math:`(A,B)` slot of :math:`N` is **structurally** zero
because the seeding term lives in :math:`M` — there is no block to
suppress, no branch to take.  :math:`M`'s ``solve`` is the block
back-substitution: System B's march first (exactly the curvilinear sweep
order), then the bulk sweep on :math:`q_A - A_{AB}\,\psi_B`.

⭐ Placement by ends is the mechanism a later partition phase
generalises to :math:`A_{ij} = R_i\,A\,J_j` for an arbitrary carrier
partition; the labelling above is its first constructor, and the
:class:`~orpheus.numerics.coupled_system.CoupledOperator` constructor's
own block type-check is what makes a mis-placed coupling
*unconstructable* rather than merely wrong.

The law, and what it is measured at
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:meth:`~orpheus.sn.splitting.Splitting.law_residual` evaluates
:math:`\bigl(M - \sum_i N_i - A\bigr)x` on a state and returns it flat,
aligned with :math:`A x`.  It is the value's **certificate against the
Problem it was minted from**: :math:`A` is read off
``system.loss``, i.e. the operator the *posing* built, never
re-assembled from the same pieces the splitting used — so the two sides
of the comparison do not share a construction.

.. list-table:: What the law reads, and why
   :header-rows: 1
   :widths: 22 26 52

   * - arm
     - residual
     - why
   * - seedless, both schedules
     - exactly ``0.0``
     - the terms are the loss's own objects and the Gauss-Seidel split
       writes **disjoint rows** — on any row carrying a boundary value
       exactly one of :math:`B_{\rm lower}`, :math:`B_{\rm upper}` is
       non-zero, so no addition is reordered.  ``[M]`` bit-exact on 6
       configurations (absorber, :math:`c = 0.5`, :math:`c = 0.9`,
       x-vacuum, at two mesh sizes) under an *exhaustive* dense
       assembly — one unit probe per degree of freedom, so every entry
       of :math:`M - N - A` is checked rather than a random subspace.
   * - carrying
     - round-off, **not** zero
     - the block-grid assembly re-associates the sums, so bit-exactness
       is not available and is not asserted.  The draw-stable statistic
       is the **relative** one: ``[M]`` :math:`\max\lvert r
       \rvert_\infty / \lVert A x \rVert_\infty = 2.087\times10^{-16}`
       (:math:`\approx 0.94\,\varepsilon`) over 40 random draws.  The
       shipped gate pins ``nulp=8`` and ``[M]`` fails at 4.

.. warning::

   ⚠ **Do not pin the carrying arm's absolute residual.**  It is a
   property of the draw, not of the splitting: ``[M]`` the same fixture
   reads :math:`3.55\times10^{-15}` at one draw and
   :math:`2.84\times10^{-14}` over forty.  Only the relative figure
   above is stable, and the gate is written against a nulp band for
   that reason.  The same caution applies to every mutation magnitude
   in this section — each is the norm of *some operator applied to some
   probe state*, so it moves with the fixture and the seed.  The gate
   docstrings carry the current numbers and re-measure them; quote
   those, not these.

The law's teeth are its mutations, and the ladder is instructive because
it shows what the law can and cannot attribute:

* **Drop every gain** (claim :math:`A = M`) — this is the historical
  σ_r-fold defect in its purest form, ERR-070, which shipped 46–56 %
  silent flux errors gated by nothing.  ``[M]`` it moves the law by
  ``1.00e-01`` seedless and ``1.83e-02`` carrying, fifteen or more
  orders above the contract.
* **Flip one gain's sign** — distinct from the drop, because the gain is
  still *present* (an arity or ``None``-block check still passes) and
  only its value moves, so this is the catcher for a convention drift.
  ``[M]`` ``3.32e-02`` seedless, ``3.66e-02`` carrying.
* **Re-introduce the ERR-056 boundary-Gauss-Seidel defect** — keep
  :math:`B_{\rm lower}` implicit *and* lag the whole :math:`B_a`.  This
  reddens the law on the Gauss-Seidel arm; the per-face partition gate
  stays green, because :math:`B_a` itself is untouched.

⛔ **The law alone cannot tell a double-label from a dropped term.**
Labelling :math:`B_{\rm lower}` *both* implicit and explicit, and
*removing* it from the implicit set, shift the residual by exactly the
same quantity — :math:`\lvert B_{\rm lower}\,x\rvert` — so the two
defects are indistinguishable by magnitude.  Separating them needs a
**structural** assertion instead: that the two label sets are disjoint
(catches the double-label) and that their multiset union is the record's
own term set (catches the drop).  That check is exact, costs nothing,
and is the only thing that discriminates.

⚠ One leg of the design has **no production home** and is a test-side
obligation: the *realization* law, that :math:`M` as derived equals the
sum of its own implicit terms.  ``law_residual`` cannot see it, because
:math:`M` is derived *from* those terms — every piece-labelling mutation
leaves it green.  Its only catcher is a perturbation of the scheduled
composite's ``apply``; without that arm the realization claim ships
unwitnessed, and the honest thing is to say so here rather than let the
law's name imply coverage it does not have.

Finally, the two labellings together are what a Mode-9 gate needs.  A
splitting must not move the fixed point, only the rate, and the honest
way to check that is to compare two *different* splittings of the same
:math:`A` on a configuration where the degeneracy that would hide an
error is broken.  ⚠ On a closed reflective diamond box :math:`A` is
exactly singular, so the two arms legitimately return different
*members* of a solution manifold and the comparison must be made on a
gauged trace or on a functional orthogonal to :math:`\ker A` — see
:ref:`sn-loss-kernel-gauge` and :ref:`sn-exit-gauge`.

The query contract: POSED operators versus USED operators
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The split above is one instance of a general rule this campaign adopted,
and it is worth stating in its own words because it decides where every
future strategy object lives:

   The Problem answers with the operators as **POSED** — original,
   unmodified, once per Problem.  The Strategy answers with the
   operators as it **USES** them — labelled, split, re-bound, lowered —
   and every Strategy value certifies itself against the Problem by an
   algebraic law.

:class:`~orpheus.sn.splitting.Splitting` is the first such value, and
the law it certifies itself by is :eq:`sn-splitting-law`.  The pattern
generalises: the boundary split certifies itself by :math:`B_{\rm
lower} + B_{\rm upper} = B_a`; a windowed gain certifies itself against
the original's own projection; a shifted factor set will certify itself
against the pencil's :math:`\mathcal A(\sigma)`.

The record the inner solve leaves behind carries **both faces in one
object**, which is why it can answer either question without a caller
guessing: :class:`~orpheus.sn.solver.InnerSolve` holds ``system`` (the
Problem's posed record), ``splitting`` (the Strategy value that ran),
``driven_gains`` (those of the value's explicit pieces that were
actually applied — the angular lifts re-bound onto the moment iterate
when the 2-D arm was windowed), and ``iterate``.  The finalize
(:ref:`sn-finalize-one-step`) reads :math:`M` and the driven gains from
it and applies the map once; because nothing re-selects a splitting the
inner already chose, the reconstruction cannot drift from the iteration
it is reconstructing.

What is deferred, and why
~~~~~~~~~~~~~~~~~~~~~~~~~~

Several things the value's shape anticipates are deliberately *not*
built.  Each is deferred for a stated structural reason, not for
effort.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - deferred
     - the structural reason
   * - a **gain grid on both arms** (one fused :math:`N` operator
       seedless as well as carrying)
     - The 2-D windowed driver re-binds the lagged angular lifts onto
       the moment domain *one by one*, and only
       :class:`~orpheus.transport.operators.angular_lift.AngularLift`
       carries an ``on_moment_domain`` re-binding.  The boundary
       operators carry none, and
       :class:`~orpheus.numerics.operator.OperatorSum` validates its
       members' domains at construction, so a fused
       ``S.on_moment_domain() + B_a`` **raises**.  Landing the grid now
       would mean a distributing ``on_moment_domain`` on
       ``OperatorSum`` and ``CoupledOperator`` plus a domain re-bind on
       two boundary classes — a new arm on four classes for one
       consumer.  Keeping the explicit part addressable *as pieces* is
       what lets the window work at all, and the grid is a derivation
       from the same primitive whenever it is wanted.
   * - a **second constructor** from a carrier partition
     - Placement by ends already generalises to
       :math:`A_{ij} = R_i A J_j`; what is missing is the partition
       object, not the mechanism.  One interface, two constructors.
   * - retiring
       :class:`~orpheus.sn.operators.scheduled_invertible.ScheduledInvertibleOperator`
     - It is the reified :math:`(L+C) - B_{\rm lower}` composite, and it
       is ruled to *retire*, not to be renamed, when the partition
       constructor lands.  This is why
       :attr:`~orpheus.sn.splitting.Splitting.implicit` is annotated by
       the **capability** it must carry — an invertible operator — and
       never by that class: the retirement then touches no signature
       here.
   * - a **policy** that chooses among splittings
     - Choosing *between* values (by spectral radius, by measured rate)
       is a different object from constructing one, and it needs at
       least two values to choose from before it can be written.  The
       transfer move it will be built on is already law-invariant
       (above).
   * - the **shift** :math:`\sigma` — an α-eigenvalue or a
       shift-invert point
     - The absorbed shift :math:`M[\Sigma_t - \sigma/v]` is a Strategy
       *lowering* of the pencil evaluated at :math:`\sigma`, not a
       property of the splitting: a shifted factor set is minted and
       then labelled exactly like this one, and its own law certifies
       it against :math:`\mathcal A(\sigma)`.
   * - the Wielandt shift :math:`-(1/k)F` as an explicit term
     - It is a term like any other once the pencil's point is chosen;
       it consumes the assignment this value carries.
   * - a schedule derived from the sweep graph's strongly-connected
       components
     - It arrives behind the same
       :class:`~orpheus.sn.loss_representation.sweep_schedule.SweepSchedule`
       call, so the labelling site does not change — the schedule is
       what the labelling *reads*, not what it *is*.

Designs that were considered and refuted
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Recording these is the point of the section: each looked reasonable, and
each fails for a reason that will still be there next time.

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - candidate
     - why it fails
   * - Store :math:`M` and :math:`N` **as fields** (the products), as
       the posed record used to
     - Then the assignment exists nowhere as data, a transfer between
       labels has no operand, and the two members are free to disagree
       with the labelling that supposedly produced them.  Deriving them
       is what makes the disagreement unspellable.
   * - Put ``Splitting`` in ``orpheus/numerics/``
     - Its one constructor's body is SN fusion — the
       ``StreamingCollisionOperator`` subtraction dispatch and the
       carrying quartet.  A numerics-layer type whose only constructor
       is method-specific is premature generalisation; the abstraction
       is earned when a second method needs it.
   * - Bind the value to the **factors alone**, not to the posed record
     - This was the adversarial review's own recommendation, and what
       shipped is wider, because two of the value's three uses need
       record members rather than leaves: the law certifies against
       ``system.loss`` (an independently-composed operator, which is
       exactly what makes the check non-tautological), and block
       placement reads ``system.space`` as its address book.  A
       factors-only value would have had to re-assemble :math:`A` from
       the same pieces it splits — comparing a construction with
       itself.
   * - Make :class:`~orpheus.sn.loss_representation.sweep_schedule.SweepSchedule`
       *the* assignment datum
     - It is an octant **order** plus a set of reflecting faces.  The
       carrying arm's assignment is not in it at all, and putting two
       rules under one name would leave the second legal carrying
       splitting (block-Jacobi, with the seeding term lagged)
       unspellable.  The schedule is an input to the labelling; the
       structural fact the labelling reads off it is just
       :attr:`~orpheus.sn.loss_representation.sweep_schedule.SweepSchedule.is_sequenced`
       — more than one octant group — with ``kind`` left diagnostic.
   * - Annotate ``Splitting.implicit`` as
       ``ScheduledInvertibleOperator``
     - "Retire, don't rename": that class is scheduled for retirement,
       and a signature naming it would have to be edited when it goes.
       The annotation names the capability instead.
   * - Discriminate the ``inner_schedule`` string where it is used
     - It was checked in two places — a repeated conditional is a
       missing type.  It is now resolved **once**, at solver
       construction, into a schedule object; the string survives on the
       entry-point signatures as sugar and nothing downstream reads it.

.. _sn-one-fission-per-problem:

ONE :math:`F` per Problem
--------------------------

The splitting section above is one half of the consumers campaign's step
2: *what the Strategy owns*.  This is the other half — *what the Problem
owns* — and it is the smaller change with the longer reach, because the
object it makes single is one of the two members of the k-eigenvalue
**pencil** :math:`(A, F)`.

A pencil is a pair of operators **on one space**.  Until 2026-09-13 the
S\ :sub:`N` forward solve and its adjoint did not have one :math:`F` to
pair with :math:`A`: they built different objects, on different spaces,
from the same data.

The defect: two mints, two spaces
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``[M]`` 2026-09-12, a counting spy on both fission factories over a
single hub (a two-region 2-group slab, ``gauss_legendre(8)``):

.. list-table::
   :header-rows: 1
   :widths: 24 40 36

   * - route
     - the object it minted
     - the space it lived on
   * - forward k-outer
     - :class:`~orpheus.transport.operators.isotropic_transfer.IsotropicFission`
       (``SNSolver.__init__``)
     - :attr:`SNProblem.bulk_space <orpheus.sn.problem.SNProblem.bulk_space>`,
       :math:`(n_g, *\text{spatial})`
   * - adjoint, **seedless** mesh
     - :class:`~orpheus.transport.operators.fission.FissionOperator`
       (``_adjoint_posing_parts``)
     - :attr:`SNProblem.full_field_space <orpheus.sn.problem.SNProblem.full_field_space>`
   * - adjoint, **carrying** mesh
     - an :class:`~orpheus.numerics.operator.OperatorProduct`
       (``_adjoint_posing_parts``, again)
     - the :class:`~orpheus.numerics.coupled_system.CoupledSpace`

Three spellings, two mints per solve, and — the part that matters —
**no object that both faces read**.  So
:math:`F_{\rm adjoint} = F^{\dagger}` was not merely unasserted: it was
*unstatable*, because the :math:`F` on the right of the sentence did not
exist anywhere the forward solve could name.  That is the same shape as
the retained-Legendre-order defect step 1 closed
(:ref:`sn-hub-retained-order`) — **a datum with three spellings is a
datum with no owner** — and it is worse here, because a spelling that
lives on a different *space* cannot even be compared to its sibling by
``==``.

⚠ Nothing in any value could see it.  Both mints read the same
:class:`~orpheus.transport.kernels.FissionKernel` pair
:math:`(\chi, \nu\Sigma_f)` and produced the same numbers; the adjoint's
:math:`k^{\dagger}` agreed with the forward :math:`k` exactly as it
should.  A duplicated *derivation* that agrees is invisible to every
value gate by construction, which is why the catcher is a **mint
count** and not a residual.

The ruling: the composite is primary
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Both objects are legitimate bindings of one datum
(:ref:`sn-fission-binding-adjoint`), so "which one is *the* :math:`F`?"
is a real question with two defensible answers.  It was ruled — open
ruling **O-1** of the step-2 design record — in favour of the
**composite**:

   :math:`F` is the angular binding
   :class:`~orpheus.transport.operators.fission.FissionOperator` on the
   full field; the scalar dyad
   :class:`~orpheus.transport.operators.isotropic_transfer.IsotropicFission`
   is a **Strategy-side reduction** of it, not a peer.

The reason is the pencil, not the fission channel.  :math:`A` is posed
on the full field (bulk :math:`\oplus` trace) because the boundary law
is a first-class sibling operator; a pencil member posed on the *bulk
alone* is not on the same space as its partner and cannot be paired with
it by any object.  The composite is therefore the member that can be
held; the scalar face is what a particular *strategy* — the k-outer's
scalar power iteration — reduces it to.

⭐ That ruling is about the long-term shape and was taken with the
re-baseline cost declared irrelevant.  What step 2 **landed** is the
bit-identical half of it (next section); moving the forward k-outer
itself onto the composite carrier is a later, principled re-baseline.

The mechanism: one cached member, two faces
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:attr:`SNProblem.fission <orpheus.sn.problem.SNProblem.fission>`
is a :func:`~functools.cached_property` on the Problem hub, minted once
through the tier-2 factory
:meth:`FissionOperator.from_solver_data
<orpheus.transport.operators.fission.FissionOperator.from_solver_data>`
with ``space=self.full_field_space``.  Everything else is a *view* of
that one object:

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - consumer
     - what it reads
     - why that face
   * - forward k-outer
       (:meth:`~orpheus.sn.solver.SNSolver.compute_fission_source`)
     - ``problem.fission.isotropic_energy``
     - the outer iterates a **scalar** flux, so it wants the rank-1
       energy dyad on :math:`(n_g, *\text{spatial})` arrays.  The face is
       a declared field of
       :class:`~orpheus.transport.operators.angular_lift.AngularLift`,
       bound at construction from the composite's own datum — a
       *theorem* of the composite, not a second binding
       (:ref:`sn-fission-binding-adjoint`)
   * - adjoint, seedless
     - ``system.factors.fission`` — *is* ``problem.fission``
     - the entry daggers it as ``F.H``; the composite is what carries
       the two metrics the Hilbert adjoint needs
   * - adjoint, carrying
     - :attr:`~orpheus.sn.coupled_system.WithinGroupSystem.production`
     - the same :math:`F`, **posed** on the coupled carrier (below)

The hub was the right home rather than the posed record because, **when
this unit landed**, the forward eigenvalue path still rebuilt that record
once per outer step.  A record-resident :math:`F` would therefore have
been re-minted every outer, and the count gate below would have read the
outer count rather than one.  ⚠ The adjoint and fixed-source paths
already built once per Problem, so the distinction was invisible on those
arms — which is exactly why a record-resident :math:`F` would have looked
correct from two of the three entry points.

⭐ **The constraint is gone, and the choice it forced was the right one
anyway.**  The campaign's next unit (C3b, the same day) made the record a
:func:`~functools.cached_property` of the hub, so it is now built ONCE per
Problem and a record-resident :math:`F` would no longer be re-minted
(:ref:`sn-the-problem-poses-its-pencil`).  The hub remains :math:`F`'s
home regardless, for a reason that does not depend on the build count:
:math:`F` is a **datum of the Problem**, not a member of the within-group
*posing* — within-group fission is zero, and the record carries the
posed :attr:`~orpheus.sn.coupled_system.WithinGroupSystem.production`
precisely as a *reference* to the hub's leaf rather than as a second
mint.  What the timing constraint bought was that the right answer was
also the only landable one.

The posed production, on either arity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:func:`~orpheus.sn.coupled_system.build_within_group_system` now poses
:math:`F` on the same carrier it poses :math:`A` on, and stores it as
:attr:`~orpheus.sn.coupled_system.WithinGroupSystem.production`.  On a
**seedless** mesh that is the :math:`1\times1` grid
:math:`[\,[F]\,]`.  On a **carrying** mesh it is a composition, not a
block matrix with holes:

.. math::
   :label: sn-posed-production-carrying

   F_{\rm posed}
   \;=\;
   \begin{bmatrix} F \\[2pt] E_F \end{bmatrix}
   \circ\;
   r_{\rm bulk} ,
   \qquad
   E_F \;=\; \mathrm{Fold}\circ F_{\text{iso}}\circ\!\int\! d\mu ,

.. (vv-status rationale) Representational identity: it states how the
   shipped builder spells the posed production on a seed-carrying mesh —
   a restriction onto the bulk member followed by a rectangular
   prolongation stack — which is a construction-site fact, not a solver
   claim.  Its verifiable content is that the posed object agrees with
   the adjoint entry's own on a random coupled state, gated by the
   ``@pytest.mark.foundation`` row
   ``tests/sn/operators/test_step2_posed_fission_anchors.py::TestLawTheHubsFissionIsWhatBothFacesRead::test_the_carrying_adjoint_poses_the_records_production``.
.. vv-status: sn-posed-production-carrying documented

where :math:`r_{\rm bulk}` is the
:class:`~orpheus.numerics.coupled_system.SystemRestrictionOperator` onto
System A (the system restriction, defined with its transpose, extension by
zero, at :ref:`coupled-block-system-restriction`) and :math:`E_F` is the
**fission ray fold** — the
kernel-generic
:class:`~orpheus.sn.operators.radial_characteristic.RadialCharacteristicEmission`
carrying ``F.isotropic_energy``.

The composition is the honest spelling of a physical fact: **fission
annihilates the ray system.**  The :math:`w = 0` closed rays carry no
quadrature weight, so nothing sources fission from them — and an
annihilated *input* column is a restriction, not a zero block.  Writing
it as a :math:`2\times2` grid with an explicit :math:`0_{BB}` was the
pre-2026-08-22 spelling, and it existed only because a
:class:`~orpheus.numerics.coupled_system.CoupledOperator` refuses an
all-``None`` column; the zero the dagger needs now falls out of the
restriction's extension-by-zero through the space's own materialization
seam (:ref:`sn-adjoint-coupled-posing`).

⭐ **This posing used to live in the adjoint entry, and now lives at the
builder.**  That is the Cardinal-Rule-2 half of the step: the forward
builder is the *only* site that knows how to put an operator on this
carrier, so a posing written anywhere else is a twin by construction —
and this one was, for as long as only the adjoint needed it.

What it cost, and what it did not
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Two routes were available for the forward k-outer, and they are not
equivalent in floating point.  Both were measured over **200 seeds**
before the ruling, on the seedless slab hub:

.. list-table::
   :header-rows: 1
   :widths: 40 18 20 22

   * - route
     - ``array_equal``
     - :math:`\max\lvert\Delta\rvert`
     - verdict
   * - read the composite's **energy face**
       (``F.isotropic_energy.apply(φ)``) — *what shipped*
     - **200 / 200**
     - ``0.000000e+00``
     - bit-identical; a pure re-homing
   * - apply the **composite** on the forward carrier
       (``F.apply(ψ)`` vs today's lift of the scalar dyad)
     - **0 / 200**
     - :math:`2.776\times10^{-17}`
       (:math:`\max` rel :math:`2.215\times10^{-16}`)
     - :math:`\le 1` nulp, draw-stable — **deferred**

The second row is not a defect; it is the *price tag* on the carrier
change O-1 rules for eventually.  The gap is an IEEE re-association and
nothing else — the lift normalises by :math:`W` before the broadcast and
the angular binding after — so it meets every criterion for a
principled-equivalent change (a named intermediate, a
structurally-independent anchor, a drift bounded by one reduction's
depth).  What makes it a *decision* rather than a formality is that the
S\ :sub:`N` regression set is already at the ULP frontier: ``[M]``
recorded at
``tests/sn/operators/test_step2_posed_fission_anchors.py``, nine of the
fourteen diamond-difference regression cases already drift 1–11 ULP
against their frozen references, so a 1-nulp shift on *every* eigen
solve would move the drift **set** rather than disappear into it.  A
re-baseline is therefore an act, not a side effect, and it is scheduled
as one.

⚠ Both figures are properties of the **fixture**, not of one draw: the
``array_equal`` claim is a 200-seed sweep, which is what licenses
pinning it at ``array_equal`` rather than at a tolerance.  A single
green reading would have licensed neither.

The gate is a COUNT, not a name
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The permanent catcher is
``tests/sn/operators/test_step2_posed_fission_anchors.py::TestRuledOneFissionPerProblem``:
it builds a forward :class:`~orpheus.sn.solver.SNSolver` **and** poses
the adjoint over one hub with both fission factories wrapped by a
counting spy, and asserts the total is ``1``.  Three properties of that
design are deliberate.

* **It counts, it does not name.**  What the hub's member is *called*
  was an open ruling when the anchor was written, so an assertion
  naming ``SNProblem.fission`` would have been a guess wearing a contract.
  The invariant is *one mint per Problem*; the attribute name is not.
* **The spy is keyed on the two ``classmethod`` factories**
  (``IsotropicFission.from_material_xs`` and
  ``FissionOperator.from_solver_data``), patched on the class so every
  module-level binding resolves through it.  This is why the hub mints
  through ``from_solver_data`` rather than calling the dataclass
  constructor directly: **a direct constructor would be invisible to the
  census**, the counter would read zero, and a gate whose instrument
  reads zero is not a gate that passed.
* **It refuses its own empty reading.**  If the counter is empty the
  helper raises rather than returning ``{}``, because an instrument that
  counted nothing carries no information about what it was pointed at —
  the positive-control discipline, spelled inside the harness.

The row is tagged ``@pytest.mark.foundation``: it is a software /
architecture invariant of the fission binding, with no theory equation
label behind it and therefore no ``verifies(...)`` marker.

Beside it sit three identity rows
(``TestLawTheHubsFissionIsWhatBothFacesRead``) that pin *which* object
each face reads — the forward's operand ``is`` the hub's energy face,
the seedless adjoint's :math:`F` ``is`` the hub's composite, and the
carrying adjoint's posing agrees with the record's ``production`` on a
random coupled state.  The count says *one was minted*; these say *the
one is the one everybody reads*.

What is deferred, and why
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 32 68

   * - deferred
     - the structural reason
   * - the k-outer's **carrier change** — the forward applying the
       composite :math:`F` on the full field
     - Ruled (O-1) as the long-term shape and measured at
       :math:`\le 1` nulp (above).  It is a principled re-baseline of
       the frozen S\ :sub:`N` references, and re-baselining is an act
       with its own three-criterion justification; bundling it into a
       re-homing would have hidden a value change inside a commit whose
       claim is that nothing moved.
   * - the :math:`1/k` **division moving off the solver**
     - :math:`1/k` is the pencil's spectral parameter, not a property
       of :math:`F` — :math:`F` is linear and the scaling is the
       *algorithm's*.  It moves when the pencil object lands and can
       carry it (``pencil.at(σ)``), not before: parking it on the
       operator now would re-create exactly the welded scaling the
       four-tier separation exists to prevent
       (:ref:`eigenvalue-posing`).

       ⚠ **Precondition discharged 2026-09-13, the row still open.**  The
       pencil landed later that day and *does* carry ``at(σ)``
       (:ref:`sn-the-problem-poses-its-pencil`) — but
       :class:`~orpheus.numerics.iteration.KEigenvalue` did not consume
       a posing yet, so the division was still on the solver.  What the
       row recorded was a *scheduled* move with a measured price, not a
       structural blocker.

       ✅ **DISCHARGED 2026-09-14** — ``KEigenvalue(posing, implicit,
       explicit, …)`` takes the question, and ``compute_keff`` is
       ``posing.rayleigh(ψ, w=1)``.  The re-baseline was the priced one:
       the loss is applied ONCE, and ``[M]`` the adjoint :math:`k` drifted
       :math:`1.22\times10^{-15}` (~5 ulp).
   * - stating :math:`F_{\rm adjoint} = F^{\dagger}` **as a theorem**
     - Step 2 made it *true* — both faces read one object — but the
       sentence is a claim about the **pencil**, and the pencil is not
       yet an object.  Once
       :class:`~orpheus.sn.coupled_system.WithinGroupSystem`'s pair is
       reified, the adjoint posing becomes ``pencil.H`` and the
       equality is a property of that type rather than a coincidence
       two entry points maintain.  Until then, the identity rows above
       pin it empirically.

       ⚠ **Half discharged 2026-09-13.**  The pair IS reified —
       :attr:`SNProblem.pencil <orpheus.sn.problem.SNProblem.pencil>`
       is an :class:`~orpheus.numerics.pencil.OperatorPencil` whose
       :attr:`~orpheus.numerics.pencil.OperatorPencil.H` daggers both
       ends — so the equality became **statable** as a property of the
       type.  It was not yet **spelled** that way: ``_adjoint_posing_parts``
       still hand-daggered the triple, because
       :class:`~orpheus.numerics.iteration.KEigenvalue` had not adopted
       the posing.

       ✅ **FULLY DISCHARGED 2026-09-17** (step 3 U2, GitHub #467).  Both
       adjoint entries now drive ``problem.eigen_posing.H()`` — ONE
       daggered posing per Problem, nullary because
       :math:`k^{\dagger} = k` needs no datum — with the seedless
       Strategy pair lifted into the :math:`1\times1` coupled grid so
       both arms share the carrier.  The equality is a property of
       :class:`~orpheus.numerics.pencil.OperatorPencil` and no entry
       spells adjoint physics.  ``[M]`` a ULP-class re-baseline:
       :math:`k_{\rm adj}` moved by rel :math:`\approx10^{-15}` against
       :math:`10^{-9}` certification gates.  The identity rows above are
       now corroboration rather than the only pin
       (:ref:`sn-adjoint-poses-a-pencil`).
   * - the record's ``production`` being **shared** between two builds
     - Two ``build_within_group_system`` calls over one hub produce two
       ``production`` objects (each composes its own restriction and
       stack) that agree by value, not by ``is``.  The leaves are
       shared — ``factors.fission`` is the hub's :math:`F` on both —
       and object identity of the *posing* arrives when the hub caches
       the record itself.

       ✅ **DISCHARGED 2026-09-13** — the hub caches the record
       (:attr:`SNProblem.system <orpheus.sn.problem.SNProblem.system>`)
       and is the builder's only caller, so every production consumer
       reads ONE ``production`` object by identity.  The literal sentence
       above stays true and stops mattering: two *direct* builder calls
       still mint two, but no production path makes them.

.. _sn-sigma-is-a-problem-datum:

:math:`\sigma_t` is a datum of the Problem
-------------------------------------------

The third act of the consumers campaign's step 2, and the one that
retires the most code.  It answers a question the two sections above
leave open: the pencil :math:`(A, F)` is a pair of operators built from
a Problem's data — so **what is the data?**  For :math:`F` the answer
was already a Problem datum (the fission kernel pair
:math:`(\chi, \nu\Sigma_f)`).  For :math:`A` it was not: the total cross
section that forms the collision diagonal :math:`C = M[\sigma_t]` could
be *rebound on a live solver*, which is a different thing entirely.

The defect: a solver you could mutate underneath its own operators
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Until 2026-09-13 the cross-section state was a solver attribute,
``SNSolver.mat_xs``, built in ``__init__``; and ``SNSolver`` carried a
method, ``rebind_cross_sections(new_sig_t)``, whose docstring named
depletion and thermal feedback as its consumers.  It wrote a new array
straight onto the field's private ``_sig_t_cell`` slot and rebuilt the
:math:`\sigma`-dependent half of the sweep cache.

Three things were wrong with it, in increasing order of reach.

#. **It mutated an object other objects already held.**  The rebind
   *re-bound* ``mat_xs._sig_t_cell`` to a new array, so an operator that
   reads through the field on every ``apply`` (:math:`S`, :math:`F`) saw
   the new :math:`\sigma_t`, while one holding a
   :class:`~orpheus.transport.fields.cross_section_field.CrossSectionField`
   built *before* the call — the collision diagonal
   :math:`C = M[\sigma_t]` — still pointed at the old array.  The
   campaign's own characterisation row measured
   this on the seedless slab anchor fixture, one seeded coupled state: a
   :class:`~orpheus.sn.coupled_system.WithinGroupSystem` built **before**
   a :math:`\times 3` rebind and applied **after** it was ``array_equal``
   to its pre-rebind self, while a freshly built one moved by
   :math:`\max|\Delta| = 3.7798`, :math:`\max\text{rel} = 0.1619`.  The
   reused system was silently stale.  The defect was unspellable only
   because
   :func:`~orpheus.sn.coupled_system.build_within_group_system` happened
   to run once per outer step.
#. **The identity question had no answer.**  The hub is the solve's save
   state (:ref:`sn-hub-retained-order`), and after a rebind it compared
   ``==`` to its pre-rebind self while representing a different problem.
   A registry keyed on the hub, a cache keyed on the hub, a
   ``same_phase_space`` pairing check — all of them would have agreed the
   two were one problem.
#. **It was guarded by an assertion that could not fire.**  ``__init__``
   carried an ``if __debug__:`` block re-deriving :math:`\sigma_t` from
   the materials and asserting it matched ``mat_xs.total_cross_section``
   — which is **false for every** :math:`\sigma_t` **that differs from
   the material derivation**.  ``[M]``
   under plain ``python``, constructing a solver over a σ-variant hub
   raises ``AssertionError: PR-INDEX-3 cell-flattening invariant broke``;
   under ``python -O`` — ORPHEUS's canonical runner — the statement is
   stripped at compile time and the same solve completes silently.  So
   the one guard on the path was simultaneously a live refusal of the
   capability C3a exists to add and a
   no-op in the suite that decides a merge: ``vv-principles`` failure
   Mode 8, and the ``coding-standards`` rule that *a bare* ``assert`` *in*
   ``orpheus/`` *is not a contract*.

The ruling (O-5 / O-6): the datum moves to the hub
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Two rulings, one shape.  **O-6** — a Problem owns **one**
:class:`~orpheus.transport.mesh.material_xs_field.MaterialXSField`:
``mat_xs`` is a :func:`~functools.cached_property` of
:class:`~orpheus.transport.mesh.material_mesh.MaterialMesh`, so every
consumer of one Problem shares one field and there is nothing to keep in
sync.  ``MaterialMesh.material_xs_field()`` — a *builder*, which minted a
fresh field on every call — retired into it, and ``SNSolver.mat_xs`` was
deleted.

**O-5** — :math:`\sigma_t` is that Problem's **datum**, not a derived
view:

.. math::

   \texttt{MaterialMesh.sigma\_t\_cell} \;\in\; \mathbb{R}^{n_g \times
   n_x \times n_y}
   \qquad\text{with}\qquad
   \texttt{sigma\_t\_cell}\Big|_{\text{default}}
   \;=\; \texttt{assemble\_cell\_xs}(\text{materials},\,
   \texttt{mat\_map}).\texttt{sig\_t}^{\mathsf T}\!\cdot\!
   \text{reshape}(n_g, *\text{spatial})

.. (Deliberately UNLABELLED.  The default derivation on the right is
   exactly :eq:`sn-cell-flatten-roundtrip`, which already carries the
   ``verifies`` witness; a second eq-label would be a twin API stating one
   law twice.  What is new here is the SLOT — that the left-hand side is a
   stored datum rather than a derived view — and that is a typing claim,
   not an equation.)

— one field, **always present**, with no ``None``-valued override and no
``is_overridden`` flag for a consumer to branch on.
:meth:`~orpheus.transport.mesh.material_mesh.MaterialMesh.with_cross_sections`
returns a **new Problem** whose datum is replaced; ``SNProblem`` re-spells it
through the same private body that serves
:meth:`~orpheus.sn.problem.SNProblem.with_scattering_order`, so
the two Problem morphisms cannot drift apart.

The datum is **normalised on the way in** — ``np.ascontiguousarray(…,
float64) + 0.0`` (which copies, and canonicalises :math:`-0.0` to
:math:`+0.0`), shape-checked against :math:`(n_g, *\text{spatial})`, and
set read-only.  That is not fastidiousness: the identity key hashes these
bytes, so re-declaring the same :math:`\sigma_t` must be the *same*
Problem, and an in-place write would desynchronise identity from content.

Identity, but not contractibility
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Where the datum lands in the two keys is the whole design:

.. list-table::
   :header-rows: 1
   :widths: 26 20 54

   * - key
     - carries σ?
     - what it answers, and what follows
   * - ``_contractibility_key``
     - **no**
     - *May two solutions' FIELDS be paired?*  The geometry, the material
       assignment, the materials' content — and, on ``SNProblem``, the
       quadrature and the scheme type.  A σ-variant leaves it untouched,
       so
       :meth:`~orpheus.transport.mesh.material_mesh.MaterialMesh.same_phase_space`
       **holds** across a depletion step: the fields pair, and
       ``Solution.compare`` / ``homogenize`` / ``condense`` keep working
       between the steps of a trajectory.
   * - ``_identity_key``
     - **yes**
     - *Is this the SAME PROBLEM?*  The contractibility key plus
       ``sigma_t_cell.tobytes()`` (and, on ``SNProblem``, the closure class
       and the clamped order).  A σ-variant compares ``!=`` and hashes
       differently, so a cache or a registry keyed on the Problem cannot
       serve one hub's table to the other by mistake.

Identity is strictly finer than contractibility
(:math:`a = b \Rightarrow a.\texttt{same\_phase\_space}(b)`), and the σ
datum is exactly the axis that separates them.

**A depletion trajectory is therefore a sequence of Problems**, not a
sequence of mutations of one.  Distinct :math:`\sigma_t^{(0)},
\sigma_t^{(1)}, \dots` give hubs :math:`P_0, P_1, \dots` with
:math:`P_i \neq P_j` for :math:`i \neq j`, and
:math:`P_i.\texttt{same\_phase\_space}(P_j)` for **every** pair
(including :math:`i = j`, and including a step that happens to leave
:math:`\sigma_t` unchanged — which is then the *same* Problem, by
construction) — the burnup step is a **morphism of Problems**, the flux of
step :math:`i` is a legitimate initial guess for step :math:`i+1` because
the fields pair, and nothing that was posed over :math:`P_i` can silently
start answering for :math:`P_{i+1}`.  That is the property the rebind
could not have: staleness is not *detected*, it is **unspellable**.

What the override reaches, and what it does not
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``mat_xs`` has four per-cell views, and after O-5 they are no longer
symmetric: ``total_cross_section`` **reads the hub's datum**, while
:math:`\sigma_a`, :math:`\nu\Sigma_f` and :math:`\chi` are still gathered
lazily from the materials through
:func:`~orpheus.data.macro_xs.cell_xs.assemble_cell_xs`.  ``[M]`` on a
:math:`\times 3` variant the first moves and the other three are
``array_equal`` to the base — which is the gate
``test_law_the_override_moves_sigma_t_and_NOTHING_else``.

A **fifth** per-cell member is neither of those two things: the
diffusion coefficient
:attr:`MaterialXSField.diffusion_coefficient
<orpheus.transport.mesh.material_xs_field.MaterialXSField.diffusion_coefficient>`
is **derived**, and since C3b-2 it is derived *from the datum* — the
subsection below is the ruling that settled it.

.. _sn-sigma-datum-diffusion-d:

The diffusion coefficient follows the datum — RULED
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. important:: **Answered.**  This subsection was a ``.. warning::``
   headed *"An OPEN question, not a settled contract: the diffusion
   coefficient does not follow"* from 2026-09-13 (unit C3a) until
   2026-09-14.  **RULED at the C3b checkpoint, fork 4 (a), and landed in
   C3b-2:** where :math:`D` is DERIVED it follows the Problem's
   :math:`\sigma_t`; where it is TABULATED it stays where it is
   tabulated.  The argument that produced the fork is preserved below,
   because it is what makes the ruling readable.

**What the open question was** (written 2026-09-13, and preserved).
:attr:`MaterialXSField.diffusion_coefficient
<orpheus.transport.mesh.material_xs_field.MaterialXSField.diffusion_coefficient>`
was a per-cell *gather* of
:attr:`Mixture.diffusion_coefficient
<orpheus.data.macro_xs.mixture.Mixture.diffusion_coefficient>`, i.e.
:math:`D = 1/(3\Sigma_{\rm tr})` with :math:`\Sigma_{\rm tr} = \Sigma_t -
\sum_{g'}\Sigma_{s,1}(g\!\to\!g')` formed from **the material's own**
:math:`\Sigma_t`.  It never read the hub's datum, so a σ-variant
:class:`~orpheus.diffusion.augmented_mesh.DiffusionMesh` posed a
diffusion problem whose **removal** term :math:`\Sigma_r = \sigma_t -
\sigma_{s,gg}` moved with the override while its **leakage** term did
not — two terms disagreeing about what :math:`\Sigma_t` is.  The two
candidate answers were (i) *declare the limit* — the override owns
:math:`\sigma_t` and nothing derived from it, the inconsistency pinned as
debt — and (ii) *close it*, which is not a one-line substitution,
because :math:`D` is **not** :math:`1/(3\sigma_t)`: the gap is the
transport correction, and "use the override" has to say what
:math:`\Sigma_{s,1}` the corrected :math:`\Sigma_{\rm tr}` is built from.

**The ruling, and why (ii) and not (i).**  :math:`D` is not an
independent datum of the problem — it is a *reading* of
:math:`\Sigma_{\rm tr}`, and :math:`\Sigma_{\rm tr}` is a reading of
:math:`\Sigma_t`.  A Problem whose :math:`\sigma_t` is stated and whose
:math:`D` is derived from a *different* :math:`\sigma_t` is not a
declared limitation, it is two Problems' data inside one hub.  The
question the fork had to answer is the one option (ii) names — *which*
:math:`\Sigma_{s,1}` — and it has a principled answer: the P1
out-scatter row sum is a property of the **material's scattering
kernel**, which the σ_t override does not touch, so it stays
per-material and only the total moves.  That splits the derivation's two
inputs cleanly, and the carve is exactly that split:

* :attr:`Mixture.p1_outflow
  <orpheus.data.macro_xs.mixture.Mixture.p1_outflow>` — new, the
  per-material :math:`\sum_{g'}\Sigma_{s,1}(g\!\to\!g')` (identically
  zero on a P0-only mixture), split out of
  :attr:`~orpheus.data.macro_xs.mixture.Mixture.transport_xs` so it can
  be consumed without the material's :math:`\Sigma_t` riding along;
* :attr:`MaterialXSField.diffusion_coefficient
  <orpheus.transport.mesh.material_xs_field.MaterialXSField.diffusion_coefficient>`
  — now **derived per cell**, :math:`D_{i,g} = 1/\bigl(3(\sigma_{t,i,g} -
  p_{1,g}[\mathrm{mat}(i)])\bigr)`, from the hub's datum
  (``mesh.sigma_t_cell``) and that per-material outflow, with a typed
  refusal when the override drives :math:`\Sigma_{\rm tr} \le 0` in any
  cell and group (*"a* :math:`\sigma_t` *override below the P1 outflow is
  not a diffusion medium"*).

``[M]`` on the reach gate's own fixture — 8 cells, 2 groups, mixture
``A``, reflective | vacuum on :math:`[0, 4]` cm, solved with
``DiffusionSolver(keff_tol=1e-10, flux_tol=1e-9)``:

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - quantity
     - base Problem
     - :math:`\times 3` σ-variant
   * - :math:`\sigma_t` (per cell, per group)
     - ``[0.5, 1.0]``
     - ``[1.5, 3.0]``
   * - :math:`p_1` outflow (per material)
     - ``[0.024, 0.045]``
     - ``[0.024, 0.045]`` — unchanged, by the ruling
   * - :math:`D = 1/(3\Sigma_{\rm tr})`
     - ``[0.70028011, 0.34904014]``
     - ``[0.22583559, 0.11280316]``
   * - :math:`1/(3\sigma_t)` (the *uncorrected* reading, for contrast)
     - ``[0.66666667, 0.33333333]``
     - ``[0.22222222, 0.11111111]``
   * - diffusion :math:`k`
     - ``0.930946184``
     - ``0.030008461``  (rel. :math:`9.678\times10^{-1}`)

Two properties of that table are the ruling's acceptance criteria, and
both are gated by
``tests/diffusion/test_sigma_variant_reach.py::test_a_sigma_variant_hub_reaches_the_diffusion_removal_term``:

#. **Bit-identical on a non-overridden hub.**  ``[M]``
   ``array_equal(mat_xs.diffusion_coefficient, gather of
   Mixture.diffusion_coefficient)`` is ``True`` — the same floats through
   the same operations, because on a hub whose datum was assembled from
   the materials :math:`\sigma_{t,i,g}` *is* the material's
   :math:`\Sigma_{t,g}`.  The re-derivation carries no arithmetic.
#. **Consistent on a σ-variant.**  Both terms now read one
   :math:`\sigma_t`.  ``[M]`` the variant's :math:`k` moves from
   ``0.029084534`` (the pre-carve, :math:`D`-frozen reading on the same
   fixture) to ``0.030008461`` — a leakage term that now shrinks with the
   override, which is the whole content of the change.

.. note:: **Two honesty notes on the numbers above, and one on the
   ruling's second clause.**

   The pre-C3b-2 text of this warning quoted ``0.26290298 →
   0.01802733`` for "an 8-cell 2-group slab, :math:`\times 3` override",
   relayed from the C3 verification delta's probe ``p11``.  ``[M]`` that
   pair does **not** reproduce on the gate's fixture — which is an 8-cell
   2-group slab with a :math:`\times 3` override — so the probe's mesh,
   mixture or boundary pair differed in a way the memo did not record.
   The rows above are re-measured here, with the fixture stated, and the
   old pair is retired rather than carried forward.  A diffusion
   eigenvalue is a property of mesh × materials × boundary; re-measure
   rather than quote.

   And the ruling has **two** clauses, of which only the first has an
   occupant.  ``[M]`` ORPHEUS ships **no tabulated** :math:`D`: every
   :math:`D` in the tree is derived from :math:`\Sigma_{\rm tr}`, and
   the legacy ``CORE1D`` ``transport`` vector is mapped onto the
   canonical ``Mixture.SigT`` / P1 moment upstream of this seam
   rather than arriving as a coefficient.  So *"stays where tabulated"*
   is a **declared future case**, not a live branch — there is no
   provenance discriminator on
   :class:`~orpheus.data.macro_xs.mixture.Mixture` and, with no tabulated
   inbound path, nothing for one to discriminate.  The clause is recorded
   so that a future tabulated-\ :math:`D` library lands as an addition
   rather than as a contradiction.

⚠ What the ruling did **not** settle is the wider question it names:
*which* of :math:`\sigma_t`'s derived quantities an override owns.
``[M]`` two further sites still read the material's :math:`\Sigma_t`
directly — the ``_gather_vector("SigT")`` calls at
``material_xs_field.py:264`` and ``:364``, inside ``MaterialXSField``'s
homogenisation and condensation projections — so a σ-variant hub's
*condensed* cross sections are still built from the materials.  That is a
separate seam with its own consumers, and it is untouched here.

The geometry table: shared by CONTENT, held by its consumers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The sweep cache is deliberately two strata (the tier boundary is derived
in :ref:`sn-curvilinear-multigroup`): a
:math:`\sigma`-free geometry table
(:class:`~orpheus.sn.sweep.cache.StreamingCoefficientCache`, Stratum 1)
and a :math:`\sigma`-bound one
(:class:`~orpheus.sn.sweep.cache.CollisionCache`, Stratum 2).  O-5 makes
that split *observable*: two σ-variant Problems are two Problems, and the
question "do they share the geometry table?" becomes answerable rather
than moot.

They do — **by construction**, because the intern
:func:`~orpheus.sn.loss_representation.geometry_cache_for` keys on the
``_contractibility_key`` (which the σ datum does not enter) crossed with
the angular-closure **class**.  The container that landed with C3a is a
small two-map object:

* a ``WeakValueDictionary`` keyed by ``(contractibility, closure class)``
  — the **sharing** half, so content-equal Problems read one table;
* a ``WeakKeyDictionary`` from hub to the tables that hub used — the
  **holding** half, so a table lives as long as a Problem that needs it.

Its ``len`` is the first map's: the number of **live distinct tables**.
So a σ-variant costs no table — ``[M]`` six σ-variants of one phase
space read ``len(_GEOM_CACHE_INTERN) == 1`` and all six get the same
object — while still taking a (weak) row in the second map, which is what
makes it a holder rather than a free rider.

⛔ **The holding half is a precondition, not an optimisation, and it was
measured.**  ``[M]`` over one 5-outer slab eigenvalue solve (8 cells,
2 groups): a weak-valued intern **with** a strong holder reads **1 build
/ 549 hits**; the same intern with **nothing** holding the value reads
**550 builds / 0 hits** — one rebuild per geometry-table resolve (then
the walk's own ``_ensure_geom_cache``, deleted at C3b-2 in favour of the
bound stratum, :ref:`sn-sigma-bound-once-at-the-operator`),
because the only reference between two sweeps was the weak one.  On that
fixture a build is ``0.169 ms`` (minimum of 15) against a ``167.8 ms``
solve (minimum of 3), i.e. **+55.4 %** wall-clock.  On a production mesh
a build costs two orders more — the P4.9b cost table in
:ref:`sn-p49b-operator-poses-with-closures` re-measures it, and 549 of
them is seconds, not milliseconds.

⚠ And **no value gate can see it**: ``[M]`` :math:`k` is bit-identical at
``0.435195214258`` under 1 build and under 550.  The instrument has to be
a **count**, and two carry it: ``test_geometry_cache_builds_exactly_once_per_mesh``
(``tests/sn/sweep/core/test_cache.py``), whose two legs pin one build
across a whole solve and across two independently posed operators over
one hub; and
``tests/sn/mesh/test_sigma_datum.py::…::test_law_the_table_is_built_ONCE_per_solve``,
which adds an explicit positive control (*the spy must observe at least
one build*) so a silently-inert spy cannot read as a clean pass.  Both
spy on ``StreamingCoefficientCache.from_mesh_and_quad``.  This is
``vv-principles`` #26 in its purest form: a route claim needs a route
instrument, because a function that redoes the work and throws it away is
indistinguishable, in its return value, from one that skipped it.

.. note:: **A count gate over the intern must start from an empty
   intern.**

   Two mechanics make a naive count unreliable, and both are properties
   of the objects rather than of the gate.  First, the intern shares by
   **content**, so a content-equal hub built by an *earlier* test can
   still be serving its table when a later one starts counting.  Second,
   a dead hub is not collected at refcount zero: the hub caches
   ``mat_xs`` and the
   :class:`~orpheus.transport.mesh.material_xs_field.MaterialXSField`
   holds ``self.mesh``, so hub and field form a **reference cycle** and
   both survive until the next cyclic collection.  The count legs
   therefore call ``_GEOM_CACHE_INTERN.clear()`` first, and the lifetime
   leg drops every reference and calls :func:`gc.collect` before asserting
   the entry is gone.

The gauge is σ-free, and that is now stated across two Problems
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:attr:`SNProblem.loss_kernel_gauge
<orpheus.sn.problem.SNProblem.loss_kernel_gauge>` is a
:func:`~functools.cached_property` over :math:`\sigma`-free data, and
the campaign pins that it stays so — otherwise a later "derive the gauge
from the pencil" simplification would rebuild it on every
:math:`\texttt{at}(\sigma)` with no value test able to notice.

Before O-5 the claim was spelled *on one solver*: rebind, then check the
hub returned the **same object**.  With :math:`\sigma` a Problem datum it
is spelled across **two hubs**, and the honest comparison is a
measurement rather than an obvious one: ``[M]``
:class:`~orpheus.sn.operators.loss_kernel_gauge.LossKernelGauge` defines
no ``__eq__`` anywhere in its MRO, so ``==`` **is object identity** —
two content-equal hubs give ``g1 == g2`` → ``False``, and a ``==`` row
would be a permanent false red.  (``apply`` is not a route either: the
gauge's domain is an
:class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace`,
which has no ``zeros()``.)  What is true, and what the re-posed gate
asserts on the seedless slab anchor, is
``np.array_equal(g1.as_matrix(), g2.as_matrix())`` → ``True``
(:math:`\max|\Delta| = 0.0`) across a :math:`\times 3` σ-variant.  The
per-hub **identity** half (one hub, one gauge object, asserted by ``is``)
stays where it was.

What retired, and what the gates are
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - retired
     - and what carries the claim now
   * - ``SNSolver.rebind_cross_sections``
     - :meth:`~orpheus.transport.mesh.material_mesh.MaterialMesh.with_cross_sections`
       — a morphism of Problems, with **no successor** on the solver by
       design.  ``test_cache.py`` #5, which used to assert the geometry
       table survived a rebind, now asserts two σ-variant hubs **share**
       it by ``is`` while each poses its own collision stratum.
   * - ``SNSolver.mat_xs``
     - the hub's ``mat_xs`` :func:`~functools.cached_property`
       (O-6).  ``TestLawOneFieldPerProblem`` counts
       ``MaterialXSField.from_mesh`` over one :math:`k`-solve and reads
       **1** — it read 2 before.
   * - ``MaterialMesh.material_xs_field()``
     - the same property.  It was a *builder*: every call minted a
       field, so "the solver's" and "the adjoint posing's" were different
       objects over one Problem.
   * - ``SNSolver.weight_norm``
     - nothing — ``[M]`` the attribute had zero readers tree-wide (the
       surviving ``weight_norm`` occurrences under ``tests/`` are
       same-named *local* variables).  The quadrature normalisation that
       matters lives in the sweep and on
       :class:`~orpheus.numerics.quadrature.Quadrature`.
   * - the ``if __debug__`` PR-INDEX-3 assert
     - ``tests/sn/mesh/test_sigma_datum.py``
       ``::TestLawTheRoundTrip`` — a real
       ``@pytest.mark.verifies("sn-cell-flatten-roundtrip")`` witness over
       both carrier tiers (:ref:`sn-cell-flattening-invariant`), plus the
       pure-storage gate promoted out of the same block at
       PR-CLEANUP-CODE §E.
   * - ``TestRecordTheStaleSigmaExposure``
     - deleted, because the exposure it characterised is now
       **unspellable**: there is no rebind for a system to go stale
       across.  A record row whose subject cannot be constructed is not a
       gate, it is a fossil.

The step's own laws land in ``tests/sn/mesh/test_sigma_datum.py``, a
``foundation`` module.  Its two *identity* classes are parametrised over
**both carrier tiers** — a bare
:class:`~orpheus.transport.mesh.material_mesh.MaterialMesh` and an
:class:`~orpheus.sn.problem.SNProblem`, because the datum is
declared at the data tier and the method mesh only extends the key: a
σ-variant is a different Problem; it shares the phase space;
re-declaring the same :math:`\sigma_t` is the *same* Problem; the stored
datum is read-only; a wrong shape is refused; and the round-trip law
above.  The field and intern classes are single-tier (SN, because they
need a quadrature and a solve): one field per Problem; a σ-variant owns
its own field; the override moves :math:`\sigma_t` and nothing else; the
operators posed over the variant move; and the intern shares across
σ-variants, is built once per solve, retains no hub in its value, dies
with its last holder, and is bounded by the number of distinct phase
spaces.

.. _sn-sigma-bound-once-at-the-operator:

:math:`\sigma` is bound ONCE, at the operator that owns it
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. important:: **Landed 2026-09-14 (C3b-2, fork 1 option (ii)).**  This
   subsection was headed *"What is still deferred"* and said that the
   :math:`\sigma`-bound stratum *"is still memoised on the hub
   (*\ ``problem._coll_cache``\ *), read back by a* ``getattr`` *with no*
   :math:`\sigma` *validation … it re-homes onto the
   StreamingCollisionOperator instance when the hub gains its posed
   record"*.  The hub gained its posed record in C3b's **first** commit
   and the memo did **not** move with it; it moves here — and not by
   re-homing the memo, which is the part the deferral got wrong.

**The defect, and why re-homing the memo would not have fixed it.**  The
1-D scan's :math:`\sigma`-bound collision table
(:class:`~orpheus.sn.sweep.cache.CollisionCache`, Stratum 2) was
memoised on the **hub** and read back with ``getattr(self.mesh,
"_coll_cache", None)`` — no :math:`\sigma` validation anywhere on the
path.  A walk handed a *second* :math:`\sigma` therefore kept marching
the **first** one's table.

``[M]`` on the witness's own fixture — a 4-cell 2-group vacuum slab at
``gauss_legendre(4)``, unit source — two sweeps at :math:`\sigma_t = 1.0`
and :math:`\sigma_t = 5.0` on ONE strategy returned ``array_equal``
results.  The statistic-free half of that reading is the one to carry:
the second answer is **identically** the first's, not merely close to it.
How *wrong* it is depends on the norm you pick, so state it: ``[M]``
against a freshly posed strategy's :math:`\sigma_t = 5` answer the
max-entrywise relative error is :math:`3.597`, the
:math:`\|\cdot\|_\infty` ratio :math:`3.288`, the :math:`L_2` ratio
:math:`2.378`.  (The verification delta's recorded ``3.573e+00`` is the
first of those three — a max-entrywise reading — which is why the other
two do not reproduce it.)  That is not a performance defect; it is a
silently wrong flux, and no value gate in the tree looked at it because
every production path happened to pose a fresh strategy per
:math:`\sigma`.

Moving that memo from the hub to the operator relocates the stash; it
does not make the stale answer **unspellable**, because the operator is
still a thing a caller can hand two :math:`\sigma`'s to.  What does is
inverting the dependency: :math:`\sigma` stops being an *argument of the
walk* and becomes a **bound stratum the walk consumes**.

**The seam.**  A walk no longer takes :math:`\sigma_t`; it takes what
:math:`\sigma` has already been bound INTO.

.. code-block:: python

   stratum = representation.bind_sigma(sig_t)        # bind ONCE
   psi, phi = representation.sweep(Q, stratum, boundary_flux)

:class:`~orpheus.sn.loss_representation.SigmaStratum` is the protocol
(one member, ``sig_t``), with exactly two realizations — one per **walk
kind**, which is the reason the protocol is not a single dataclass:

.. list-table::
   :header-rows: 1
   :widths: 26 30 44

   * - realization
     - carries
     - who binds it, and why that shape
   * - :class:`~orpheus.sn.loss_representation.RawSigmaStratum`
     - :math:`\sigma_t`
     - the base ``_LossRepresentation.bind_sigma``, and
       :meth:`ScanMarch.bind_sigma
       <orpheus.sn.loss_representation.ScanMarch.bind_sigma>` on a
       multi-D mesh.  The wavefront walks read :math:`\sigma` **inside
       the cell update**, so there is nothing to precompute
   * - :class:`~orpheus.sn.loss_representation.ScanStratum`
     - ``(geom, coll, sig_t)``
     - :meth:`CumprodScan.bind_sigma
       <orpheus.sn.loss_representation.CumprodScan.bind_sigma>` and
       ``ScanMarch.bind_sigma``'s 1-D branch.  The Blelloch closed form
       marches a **precomputed** chain, so binding σ means building
       Stratum 2 against the interned Stratum 1

Three consequences follow, and each is a deletion rather than an
addition:

#. **The walk's own ensure-path is gone.**  ``_OneDimScanWalk`` reads
   ``stratum.geom`` / ``stratum.coll`` directly; its
   ``_ensure_geom_cache`` and ``_ensure_coll_cache`` methods are
   **deleted**, and with them the ``getattr`` that made the staleness
   spellable.
#. **Handing a raw stratum to the scan is a typed refusal**, not a
   miscompute: the private ``_scan_stratum`` narrow raises
   ``TypeError: the 1-D scan needs a ScanStratum (geometry + collision
   tables bound for one σ) — bind σ through the scan strategy's
   bind_sigma``.  A wavefront stratum has no tables to scan with, and
   the type says so.
#. **The solver's cache block is retired.**  ``SNSolver.geom_cache`` and
   ``SNSolver.coll_cache`` are **deleted**, together with the
   ``problem._coll_cache = …`` stash in ``SNSolver.__init__``.

**Who holds the geometry table now.**  The operator does, through the
stratum:
:attr:`StreamingCollisionOperator.sigma_stratum
<orpheus.sn.operators.streaming.StreamingCollisionOperator.sigma_stratum>`
is a :func:`~functools.cached_property` returning
``self.loss_representation.bind_sigma(self.sigma)`` — so on a 1-D mesh
the operator holds a ``ScanStratum``, which holds the interned Stratum-1
table strongly.  Every ``solve`` / ``solve_transpose`` on that operator
passes the SAME stratum, which is the count claim of
:ref:`the geometry-table section <sn-sigma-is-a-problem-datum>` restated
at its new owner: one table per Problem, one binding per operator.

.. note:: **The sole-guarantor warning that stood here is discharged.**

   Between C3a and this commit, ``SNSolver.geom_cache`` looked like dead
   state — ``[M]`` read at exactly one place, three lines below its own
   assignment — and this page carried a ``.. warning::`` saying *"do not
   retire it as dead"*, because it had been silently promoted to the
   **only strong holder** of the weak-valued intern's table, and dropping
   it costs the 550-builds-per-solve regression measured above with every
   value gate green (``coding-standards``' sole-guarantor shape).  The
   warning named the event that would retire it — *"it stops being the
   tree's only holder when the σ stratum re-homes onto the operator"* —
   and that event is this commit.  The slot is now deleted; the holder is
   ``sigma_stratum``.

The §6c witness lands with the seam:
``tests/sn/sweep/core/test_cache.py::test_two_sigmas_on_one_strategy_give_two_answers``
sweeps one strategy at two :math:`\sigma` and asserts (a) the two answers
**differ** and (b) the second equals a freshly posed strategy's, ``[M]``
bit-identically.  Leg (a) is the one that was red pre-carve; leg (b) is
what makes it a *correctness* gate rather than a "something changed"
gate.  ``test_cache.py``'s intern row additionally asserts ``not
hasattr(hub, "_coll_cache")`` — the retired memo cannot come back by
accident.

⚠ **One sibling memo survives, and it is not the same object.**
``_pole_mirror_cache`` — the :math:`r = 0` coupled-pole mirror pairing —
is still a hub attribute.  It is :math:`\sigma`-**free** (it is derived
from the quadrature's mirror motion), so it carries none of the
staleness this seam removes, and it is out of this unit's scope.  The
route gate's ``_MEMO_SLOTS`` tuple in
``tests/sn/operators/test_operator_feeds_the_walk.py`` now lists it and
``_geom_cache`` alone.

.. _sn-the-problem-poses-its-pencil:

The Problem poses its pencil, and the solver holds nothing
-----------------------------------------------------------

The fourth act of the consumers campaign's step 2, and the one the other
three were clearing the way for.  :math:`F` became a Problem datum
(:ref:`sn-one-fission-per-problem`); :math:`\sigma_t` became a Problem
datum (:ref:`sn-sigma-is-a-problem-datum`); the splitting left the posed
record for a Strategy value
(:ref:`sn-splitting-is-a-strategy-value`).  What was left was the step
the whole chain exists for: **a Problem's last step is its pencil**, and
nothing before this one could hold it, because a pencil is a pair
:math:`(A, M)` on ONE space and S\ :sub:`N` had no single object to pair
with :math:`A`.

The general theory of the object — what a source-driven and an
eigenvalue system each ARE as pencils, the four :math:`(M?, q?)` cells,
the degree contract, why the shape with ``Optional`` fields lost — is
:ref:`the-operator-pencil` on the operator-algebra page, because the type
is method-agnostic and the infinite-medium solver poses the same one.
This section is the S\ :sub:`N` half: what the hub now carries, what the
solver stopped carrying, and what the change did and did not move.

The chain, and where it stops
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   hub          = SNProblem.from_material_mesh(...)       # the Problem
   record       = hub.system                           # posed: space, factors, loss, production
   pencil       = hub.pencil                           # OperatorPencil(record.loss, record.production)
   question     = hub.eigen_posing                     # EigenPosing(pencil, K_MAP)

Each of the three is a :func:`~functools.cached_property` over immutable
data, so ``hub.system is hub.system`` and the pencil is over the hub's
OWN record, not a fresh build: ``[M]``
``hub.pencil.lhs is hub.system.loss`` and ``hub.pencil.rhs is
hub.system.production``, both by identity.

The chain stops at the question.  It does **not** continue into a
resolvent, an inverse, a splitting, a schedule or a tolerance — and that
is gated rather than merely stated.
``tests/sn/architecture/test_posing.py`` walks every callable on
``SNProblem(...) → .system → .pencil → .eigen_posing`` and asserts that none
of them accepts a Strategy token (``inner_solver``, ``inner_schedule``,
``max_iter``, ``max_inner``, ``tol``, ``inner_tol``, ``restart``,
``corrector``, ``preconditioner``, ``n_dof``, ``initial_guess``).  The
gate is non-tautological in the way that matters: it holds the chain as
an explicit LIST of nine callables, so a renamed member empties the loop
loudly instead of silently, and the mutation that reds it is adding
``inner_schedule: str = "jacobi"`` to
:func:`~orpheus.sn.coupled_system.build_within_group_system`.

Its sibling row asserts the other half of what "determined by the
generating data" means: two content-equal hubs pose **equal** records,
and each hub's pencil is its own.

What the seam was, and what deleting it bought
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Until this unit, :func:`~orpheus.sn.coupled_system.build_within_group_system`
took two keyword arguments, ``scattering_op=`` and ``n2n_op=``, whose
purpose was to let the *solver* hand its own cached :math:`S` and
:math:`N_{2n}` into the build so that the two copies would not be minted
twice.  That is a caching seam wearing an injection interface, and it has
a consequence nobody wanted to state: while it existed, the Problem's
:math:`A` was **not a function of the Problem's data** — it was a
function of the data *and of whatever the caller passed*.

``[M]`` the census, with its predicate stated because the number depends
on it: **34** keyword injections across **14** files counting *consumer*
call sites (AST ``keyword`` nodes named ``scattering_op`` / ``n2n_op`` on
any ``Call``, over the tracked tree, excluding the solver's own 8
forwarding calls); **43** across **16** counting every such node,
untracked probes included.  The interesting split is inside the 34: **33**
passed back the hub's own operator — the seam was paid for and not used —
and exactly **one** minted a foreign P0 :math:`S`, a reference
construction inside ``tests/sn/operators/test_psi_half_coupling.py``,
which now poses on ``sn.system`` like every other consumer.  Both keywords
are deleted, together with ``build_coupled_system``'s long-dead
``scattering_order`` parameter.

The solver's two matching slots went with them.  ``SNSolver.scattering_op``
and ``SNSolver.n2n_op`` are **deleted**, and every read migrated to
``problem.system.factors.scattering`` / ``.n2n``.  Two censuses, because
the predicate matters: ``[M]`` **199** reads across **31** files under the
carve's own predicate (AST ``Attribute`` nodes whose *receiver* is
``SNSolver``-derived), and ``[M]`` **216** across **33** for every
``Attribute`` node of either name whatever the receiver — the extra 17
being locals, parameters and the injection keywords that retired with the
builder's signature.  With ``fission_op`` already gone at unit C2, the
solver now caches **no operator at all** — the head of this chapter says
what each slot became.

The build count, and the anchor that flipped
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The measurable claim of this unit is a **count**, not a value.

Before it, the forward eigenvalue path called the builder inside the
outer loop, so it re-posed the whole within-group system — the streaming
composite, the scattering leaf, the boundary operator, the fission
production — on **every outer step**.  The adjoint and fixed-source
paths already built once per solve, which is exactly why the defect was
invisible from two of the three entry points: any gate written on those
arms reads 1 either way.

``[M]``, from the campaign's pre-carve anchor (a counting spy bound in
both modules that call the builder, on a fissile slab that takes more
than one outer step so that "once per outer" and "once per Problem" are
different numbers):

.. list-table:: Builder calls per solve — the pre-carve anchor's own rows
   :header-rows: 1
   :widths: 40 14 23 23

   * - solve
     - ``n_outer``
     - builds, before
     - builds, after
   * - 1-D slab eigenvalue, source iteration
     - 4
     - 4
     - **1**
   * - 1-D slab eigenvalue, Krylov
     - 4
     - 4
     - **1**
   * - 1-D sphere eigenvalue (carrying), SI
     - 5
     - 5
     - **1**
   * - 2-D Cartesian eigenvalue, SI
     - 3
     - 3
     - **1**
   * - 1-D slab fixed source
     - n/a
     - 1
     - 1
   * - 1-D slab adjoint eigenvalue
     - 4
     - 1
     - 1

The "before" column is a **fixture** reading (those ``n_outer`` values
belong to the anchor's specific problems); the "after" column is a
**law**, and the row that asserts it re-measures both halves itself —
``tests/sn/architecture/test_step2_terminal_object_anchors.py``, whose
docstring carries the table above and is the thing to re-run rather than
this page.

The anchor that recorded the *pre-carve* behaviour asserted ``spy.calls
== n_outer`` — deliberately the MECHANISM rather than a fixture number,
so that a naive "move the call into a method" would still read
``n_outer`` and still red.  It was deleted with the carve, and its ruled
successor
(``TestRuledTheBuildIsOncePerProblem.test_ruled_eigen_builds_once``) lost
its ``xfail(strict=True)`` marker in the same commit — the marker's own
reason named the event that would retire it, and the event happened.

⭐ **No value moved.**  The forward :math:`k` path is bit-identical
across the change: ``[M]`` the anchors' slab reads
:math:`k = 0.4351952142580926` before and after, and the escalated
regression set is unchanged.  That is the whole design intent of building
the same operators once instead of :math:`n_{\rm outer}` times — the
operators were already functions of Problem data, so re-minting them was
pure waste, not a different answer.

.. _sn-the-pencil-reaches-production:

What the unit's SECOND commit landed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. important:: **Answered.**  This subsection was headed *"What is NOT in
   this commit"* and listed four pieces of the design as *written and not
   shipped*, so a reader would not go looking for them.  All four landed
   on 2026-09-14 in the unit's second commit (C3b-2); the list below is
   the same four, in the same order, with what each one turned out to
   be.  The fourth — the :math:`\sigma`-bound sweep memo — is
   :ref:`sn-sigma-bound-once-at-the-operator`, and is not repeated here.

**1 —** :class:`~orpheus.numerics.iteration.KEigenvalue` **consumes the
posing.**  The signature is now ``KEigenvalue(posing, implicit,
explicit, …)``: the *question* and the *Strategy's* two operators, where
it used to take an operator triple ``(A, S, F)`` and re-derive the
question from the argument order.  The three estimators read the pencil
the posing carries —
:meth:`~orpheus.numerics.iteration.KEigenvalue.compute_fission_source`
and
:meth:`~orpheus.numerics.iteration.KEigenvalue.compute_production_rate`
read ``posing.pencil.rhs``, and
:meth:`~orpheus.numerics.iteration.KEigenvalue.compute_keff` is
``posing.rayleigh(ψ, w=1)``, the one Rayleigh body of
:ref:`the-operator-pencil`.

That re-signature is the ULP-level re-baseline the deferral was written
for, and it lands with its measurement:

.. math::

   k \;=\; \frac{\sum (F\psi)}{\sum\bigl((A-S)\psi\bigr)}
   \qquad\text{replaces}\qquad
   k \;=\; \frac{\sum (F\psi)}{\sum (A\psi) - \sum (S\psi)} .

The **only** change is that the loss is applied **once**.  Mathematically
they are the same number; in IEEE-754 they are a re-association, and
``[M]`` 1 of 40 random draws is bit-identical (the C3 verification
delta's §D.3 reading, on its own fixture).  The rate is not a universal
— it is governed by how much cancellation the subtraction does, i.e. by
the scattering ratio :math:`c = \sum(S\psi)/\sum(A\psi)`.  ``[M]`` an
independent numpy re-derivation over 400 draws per row (:math:`n = 64`,
dense random :math:`A`, :math:`S = c\,A` jittered :math:`\pm 10\,\%`):

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - :math:`c \approx \sum(S\psi)/\sum(A\psi)`
     - draws where the two spellings are bit-identical
   * - :math:`0` (the ``ZeroOperator`` posture)
     - **400 / 400** — structurally, there is nothing to re-associate
   * - :math:`0.5`
     - 145 / 400
   * - :math:`0.9`
     - 43 / 400
   * - :math:`0.99`
     - **3 / 400**

So "1 of 40" is exactly where a scattering-dominated transport fixture
belongs on that curve, and the ``S = ZeroOperator`` row of
``tests/numerics/test_estimators_as_functionals.py`` is bit-identical
**by construction**, not by luck.

``[M]`` the end-to-end drift on the path the deferral named — the
adjoint :math:`k`.  On the anchors' two-region 2-group slab
(``gauss_legendre(8)``, 4 + 4 cells, ``keff_tol = 1e-10``,
``inner_tol = 1e-11``), driving :func:`~orpheus.sn.solver.solve_sn_adjoint`
with the new estimator against the retired association *in the same
process*:

.. list-table::
   :header-rows: 1
   :widths: 46 54

   * - reading
     - value
   * - forward :math:`k`
     - ``0.4351952113470525``
   * - adjoint :math:`k`, new — :math:`\sum(F\psi)/\sum((A-S)\psi)`
     - ``0.4351952113428993``
   * - adjoint :math:`k`, retired — :math:`\sum(F\psi)/(\sum(A\psi)-\sum(S\psi))`
     - ``0.4351952113428981``
   * - drift
     - :math:`1.22\times10^{-15}` absolute, :math:`2.81\times10^{-15}`
       relative — about 5 ulp
   * - :math:`|k^\dagger - k|`
     - ``4.153e-12`` new against ``4.154e-12`` retired — the
       forward/adjoint agreement is **unmoved**

That last row is the one that matters: the re-baseline is four orders
below the convergence residual it sits inside, so every adjoint
certification row stays green under its own tolerance, and the
:math:`k^\dagger = k` identity is not degraded.  The adjoint chapter
carries the posing's new spelling (:ref:`sn-adjoint-daggered-posing`).

⭐ **Two numerical facts the single-sourcing exposed**, both of which
had to be *designed in* rather than discovered afterwards, because
either one silently breaks the bit-identity the gate asserts:

* the pairing must reduce with :func:`numpy.sum`, not a BLAS dot.
  ``[M]`` on a 4096-entry vector ``np.vdot(1, y)`` and ``np.sum(1 * y)``
  differ in the last digit (``2047.5331122904424`` vs
  ``2047.5331122904417``) — pairwise summation is not the dot product's
  accumulation order, and the method-tier estimators take the coordinate
  sum.  :func:`orpheus.numerics.posing` therefore spells the pairing
  ``float(np.sum(w * y))``;
* the spectral map must state :math:`\lambda` as **one** division.
  ``[M]`` ``1/(41/6)`` is ``0.14634146341463417`` and ``6/41`` is
  ``0.14634146341463414`` — two divisions are not one.  Hence
  :attr:`SpectralMap.of_quotient
  <orpheus.numerics.posing.SpectralMap.of_quotient>` beside
  ``forward``: ``K_MAP.of_quotient(a, m) = m / a`` is the ratio
  *directly*, where ``forward(a / m)`` would round twice.  (``41`` and
  ``6`` are that gate's own :math:`\sum(A\psi)` and :math:`\sum(F\psi)`
  — the fact is not abstract, it is the row.)

**2 —** :meth:`SNProblem.source_posing(q)
<orpheus.sn.problem.SNProblem.source_posing>` **is a member.**
It returns ``SourcePosing(self.pencil.at(1.0), q)`` — **always**
``at(1)``, with no discrimination on the datum (RULED 2026-09-13, fork 2
(a)).  On a non-fissile hub the production is the zero dyad, so the
member *is* the pure transport operator in value; a branch on "does this
hub carry :math:`\nu\Sigma_f`?" would be a second spelling of a fact the
algebra already states.

One thing the implementation had to learn, and it is worth carrying: the
posing's **ends law** refused the first attempt.  ``[M]`` at the probe, a
bare fixed-source right-hand side is a ``TimedFullField`` on the
``FullFieldSpace`` while ``pencil.at(1)``'s ends are the ONE-system
``CoupledSpace`` — equal in content, different as spaces — so
:class:`~orpheus.numerics.posing.SourcePosing`'s ``__post_init__``
rejected the pair.  That is the guard working, not a nuisance: the fix is
to **lift**, not to loosen.  ``source_posing`` wraps a bare member as
``CoupledField(systems=(source,))``; a rhs that is already coupled passes
through.

**3 — the subcritical multiplying source ships**, as
:func:`~orpheus.sn.solver.solve_sn_multiplying_source`, and it is the
:math:`(M, q)` cell's production witness.  It has its own subsection
below, because it is a capability and not only a gate.

.. _sn-subcritical-multiplying-source:

The subcritical multiplying source — the :math:`(M, q)` cell, in production
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A **fixed source in a multiplying medium** — an accelerator-driven
system, a subcritical start-up source, a source-driven experiment — is
the fourth cell of the :math:`(M?, q?)` table
(:ref:`the-operator-pencil`), and until this commit that cell was argued
and not shipped.  The physics is

.. math::
   :label: sn-multiplying-source

   \bigl(L + C - S - N_{2n} - B - F\bigr)\,\psi \;=\; q ,

.. (vv-status rationale) The equation the (M, q) cell poses — the
   within-group loss MINUS the fission production, driven by an external
   source.  A posing identity rather than a solver claim: it states which
   operator the composition ``SourcePosing(pencil.at(1), q)`` builds.  Its
   verifiable content is the composition law row
   tests/sn/solve/test_subcritical_multiplying_source.py::
   test_the_composition_is_the_loss_minus_the_production (bit-identical
   against ``loss.apply(x) − production.apply(x)`` on a seeded coupled
   state, with a non-triviality positive control) together with the 0-D
   closed-form row; the physical content of a multiplying solve is
   anchored by the closed form, not here.
.. vv-status: sn-multiplying-source documented

i.e. exactly :math:`\mathcal{A}(1)\,\psi = q` — the pencil evaluated at
the physical :math:`\sigma = 1`, handed to a source posing.  In code that
is one line, and it is the hub's:

.. code-block:: python

   posing = problem.source_posing(q)        # SourcePosing(pencil.at(1), q)

**Admissibility is spectral, and the DRIVER certifies it.**  The problem
is well posed iff the medium is subcritical,

.. math::

   \rho\bigl(A^{-1}F\bigr) \;=\; k_{\rm eff} \;<\; 1 ,

because :math:`(A - F)^{-1} = A^{-1}\sum_{n\ge0}(FA^{-1})^n` converges
exactly on that condition — the Neumann series *is* the fission-chain
sum, and at :math:`k \ge 1` the chain does not terminate and no positive
steady solution exists.  That is a property of the resolvent, so nothing
on the Problem side can answer it: a posing constructs freely, and
:func:`~orpheus.sn.solver.solve_sn_multiplying_source` runs the hub's own
k-solve **first** and refuses with a typed
:class:`~orpheus.sn.solver.SupercriticalSourceProblem` naming the
measured :math:`k`.  It costs one extra eigen solve per multiplying
solve, and it is **exact** rather than a bound.

⚠ **The predicate is** :math:`\rho(A^{-1}F) < 1`, **not positive-stability
of** :math:`A - F`.  Those coincide for an M-matrix and they do not
coincide here: ``[M]`` positive-stability never holds on the S\ :sub:`N`
composite, because the trace block contributes eigenvalues :math:`-1`
(the same structure :ref:`sn-loss-kernel-gauge` analyses).  A guard
written on the "obvious" criterion would refuse every well-posed
multiplying problem in the tree.

**The lowering: the production is ONE MORE explicit gain, lagged.**  The
Strategy does not invert :math:`A - F`.  It keeps the same splitting the
pure fixed-source path uses and adds :math:`F` to the lagged side:

.. math::

   M\,\psi_{n+1} \;=\; N\,\psi_n \;+\; F\,\psi_n \;+\; q ,

which is
:func:`~orpheus.sn.solver._within_group_si`'s ``extra_gains=`` channel —
the same door the eigenvalue finalize's gain list uses, so no new
iteration body exists.  ``[M]`` the gain is posed on the **arm's own
carrier**: ``system.production`` on a carrying (curvilinear, seed-bearing)
mesh, ``system.factors.fission`` on a seedless one.  Source iteration
only: the Krylov arm's preconditioner is the *pure transport* resolvent,
and composing it with the fission lag is a later step.

Because the entry *is* that lowering, it inherits the fixed-source path's
exits wholesale — the same :ref:`exit gauge <sn-loss-kernel-gauge>`, the
same convergence certificate, the same
:class:`~orpheus.sn.solution.Solution` contract.
``tests/sn/solve/test_every_entry_gauges_its_trace.py``'s entry ledger
records it as *not separately exercised* for exactly that reason, naming
the shared path the fixed-source rows already cover — a **declared
inheritance**, not a coverage gap.

.. warning::

   ⛔ **That inheritance was NOT wholesale, and the gap was silence.**
   Until 2026-09-17 this entry returned the private SI arm's
   :class:`~orpheus.sn.solution.Solution` **directly**, so it never
   reached the two warnings its four siblings emit — both of which are
   deliberately hoisted into the PUBLIC entry so ``stacklevel=3`` blames
   the caller (#340 N4.7).  A truncated multiplying solve and a
   gauge-singular one both came back silent while the same hub solved
   through :func:`~orpheus.sn.solver.solve_sn_fixed_source` warned.  The
   defect, how it was *gated as correct*, and its catcher are
   :doc:`ERR-086 </theory/verification/error_catalog>`; the repair is one
   more thing this entry now inherits, and the emission-site count that
   pinned the old inventory reads **8** rather than 7.

   ⭐ And the entry now records what it measured to admit itself: the
   admissibility :math:`k` rides the Solution's certificate as
   ``Certified(k, "the hub's k-solve at keff_tol=…")`` — **with the
   tolerance it was measured at**, because a bare number would be a
   measurement without its configuration.  Until step 3 the driver
   measured that :math:`k`, used it, and threw it away, so a consumer
   holding the Solution could not tell a certified-subcritical answer
   from an unchecked one.

**The exit certificate is posed on the equation that was SOLVED.**  This
is the subtlety the lag creates.  The driver's certificate evaluates
:math:`r = A\psi - q_{\rm certified}`, and a lagged gain is part of the
operator, so the gain must re-enter the certified right-hand side at the
converged iterate:

.. math::

   q_{\rm certified} \;=\; q \;+\; \sum_i G_i\,\psi_{\rm conv}
   \quad\Longrightarrow\quad
   r \;=\; A\psi - \Bigl(q + \sum_i G_i\psi\Bigr)
     \;=\; \Bigl(A - \sum_i G_i\Bigr)\psi - q ,

which is :eq:`sn-multiplying-source`'s residual.  Certifying against the
bare :math:`q` would have measured the *pure transport* residual of a
multiplying solution and raised
:class:`~orpheus.sn.solver.ConvergenceCertificateError` on
every converged run.

**Convergence rate, and a budget that is not a refusal.**  The lagged
iteration's contraction is governed by :math:`k`, so the amplification
:math:`1/(1-k)` is also the iteration-count scale.  ``[M]`` on the
witness's near-critical fixture (:math:`L = 4` cm, :math:`k =
0.907457573`, so :math:`1/(1-k) = 10.81`) the solve needs **4300** inner
iterations at ``inner_tol = 1e-12``, where the budget derived from the
tolerance alone is **1961**.

.. warning:: A starved inner is a **budget**, not a refusal.  The default
   :func:`~orpheus.numerics.convergence.default_iteration_budget` is
   derived from the *tolerance*, and it does not know :math:`k`; near
   criticality it under-counts by the dominance factor.  The witness
   passes ``max_inner=6000`` on that row deliberately, and a user solving
   a near-critical driven system should do the same rather than read the
   truncation as a physics refusal (#340's convergence contract makes it
   audible either way — the exit warns and
   :meth:`Solution.converged <orpheus.sn.solution.Solution.converged>`
   reads ``False``).

**The witness** is ``tests/sn/solve/test_subcritical_multiplying_source.py``
(``l1``), six rows over a two-region 2-group slab
(mixture ``A`` | mixture ``B``, 4 + 4 cells, ``gauss_legendre(8)``):

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - leg
     - what it pins
   * - convergence ×2
     - :math:`L = 2` (:math:`k = 0.435195214`) and :math:`L = 4`
       (:math:`k = 0.907457573`) both converge; the production exit
       CERTIFICATE (above) is what asserts the balance closes, and it
       **raises** rather than reporting
   * - cone monotonicity
     - :math:`\psi_{\rm mult} > \psi_{\rm pure}` cell-wise against
       :func:`~orpheus.sn.solver.solve_sn_fixed_source` on the same
       :math:`q` — fission adds neutrons (Krein–Rutman / the Neumann
       series' positivity), and the peak ratio exceeds 3 on the
       near-critical row
   * - the refusal
     - :math:`L = 8` reflective|reflective, ``[M]`` :math:`k =
       1.374233987`, raises
       :class:`~orpheus.sn.solver.SupercriticalSourceProblem` with the
       measured :math:`k` in the message
   * - the composition law
     - ``SourcePosing(pencil.at(1), q).operator.apply(x)`` equals
       ``loss.apply(x) − production.apply(x)`` **bit-identically** on a
       randomly seeded coupled state, with a positive control asserting
       the state is non-trivial
   * - the 0-D closed form
     - the L0 oracle, below

``[M]`` **the 0-D closed form.**  On the infinite-medium pencil built
from mixture ``A`` with :math:`\nu\Sigma_f` and :math:`\Sigma_f` scaled
to make the medium subcritical, :math:`(A - F)^{-1}\mathbf{1}` has an
exact rational value and the sign flip at criticality is the refusal
leg's 0-D control:

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

   * - scale
     - :math:`k_\infty`
     - :math:`(A-F)^{-1}\mathbf{1}`
   * - :math:`0.4`
     - ``0.75``
     - ``[60, 70]`` — strictly positive, the physical solution
   * - :math:`0.6`
     - ``1.125``
     - ``[-146.667, -136.667]`` — **negative**: the operator is still
       invertible, and its solution is not a flux

The second row is the whole argument for the driver's refusal in
miniature.  :math:`A - F` does not become *singular* the moment
:math:`k` crosses 1 — it becomes **indefinite**, so a linear solve
happily returns a vector, and that vector is not a neutron flux.  A guard
keyed on invertibility would pass; only the spectral predicate refuses.

.. note:: :func:`~orpheus.sn.solver.solve_sn_fixed_source` is
   **unchanged**, and step 3 (2026-09-17) took the decision that was open
   here rather than leaving it open: on a fissile hub it still poses the
   *pure transport* operator, so every existing caller is bit-identical,
   and the multiplying question stays reachable only through **its own
   entry** — never through a boolean flag on the other one
   (``coding-elegance``: a flag parameter that selects between two
   operators is a missing type, and here the type already exists — it is
   the posing).

   What step 3 added is that the choice is now **recorded**.  The plain
   entry poses ``SourcePosing(hub.pencil.at(0.0), q)`` itself, on the
   Strategy side, because suppressing fission on a fissile deck is that
   entry's *modelling choice* and not a datum of the generating data
   (RULED 2026-09-14); the hub's own
   :meth:`~orpheus.sn.problem.SNProblem.source_posing` keeps
   naming the physical member at :math:`\sigma = 1`.  The suppression
   costs nothing to spell — ``[M]`` ``pencil.at(0.0) is pencil.lhs is
   system.loss`` by **object identity**, so "was fission suppressed?" is
   answered by the recorded posing without a second operator existing
   anywhere.  Two Solutions over one hub are therefore distinguishable by
   the question they answer, which they were not before
   (:ref:`sn-solution-carries-its-posing`).


.. _sn-finalize-one-step:

The returned angular flux — one step of the map the iteration drove
-------------------------------------------------------------------

:class:`~orpheus.sn.solution.Solution` carries **one** state, and the flux
members are readings of it — but that state is not what the power
iteration exchanged, so it has to be *reconstructed*, for two independent
reasons:

* the outer iteration's contract is scalar — ``[M]`` all **five** members
  of :class:`~orpheus.numerics.eigenvalue.EigenvalueSolver` exchange only
  the method's ``Carrier``, which for S\ :sub:`N` is the bare
  ``(n_g, *spatial)`` scalar flux, so no per-ordinate field crosses the
  boundary at all; and
* on the 2-D Cartesian windowed arm the within-group iterate is not
  per-ordinate at all — it is the harmonic-moment composite
  (:ref:`sn-angular-windowing-honest-scope`), and the user-facing
  :math:`(N, n_g, n_x, n_y)` field has to be built from it.

⛔ **Until 2026-09-17 this section opened "**\ ``Solution`` **ships TWO
flux members, and only one of them is what the power iteration
converged"** — ``scalar_flux`` was the outer's converged :math:`\phi`,
stored, and ``angular_flux`` the reconstruction beside it.  Two
representations of one quantity, each with its own space check.  Since
the consumers campaign's step 3 the reconstructed state is stored WHOLE
and :attr:`~orpheus.sn.solution.SolutionBase.scalar_flux` is
:math:`\int\psi\,d\Omega` of its cell-average moment
(:ref:`sn-solution-carries-its-posing`), so the pair cannot disagree by
construction.  What the derived reading costs is measurable and small:
``[M]`` 2026-09-17, against the #448 artefacts as frozen BEFORE the step,
``max rel`` :math:`2.19\times10^{-11}` (``cart2d_L0``) and
:math:`2.94\times10^{-11}` (``slab_vac_L0``), with :math:`k`
**bit-identical** (:math:`\lvert\Delta k\rvert = 0`) on both; the pins'
hard band is :math:`10^{-8}`, i.e. 457× and 340× of headroom.  The 16
finalize ``.npy`` and the 10 eigen regression ``.npz`` were then
**re-baselined with their generators** on the step's own tree (U2e, the
commit after the carve): a drift the bit-identity instrument
(``-W error::…DriftWarning``) would carry forever masks any LATER drift on
the same cases, so the spine was returned to a zero reading (``[M]`` 0 of 117
rows drifting after the re-baseline; the per-case table before it is in the
U2e commit message and the campaign memo's §9.3).
That gap is the same quantity as the gauge displacement measured two
sections below — the polished :math:`\psi` sits off the power
iteration's own scalar by the outer residual — so it shrinks with the
tolerances rather than being a floor.

**The reconstruction is one application of the splitting map — not a
solve.**  A within-group splitting writes the loss operator as
:math:`A = M - N` and iterates

.. math::
   :label: sn-finalize-map

   \psi^{(j+1)} \;=\; G\bigl(\psi^{(j)}\bigr)
   \;=\; M^{-1}\Bigl(q + \sum_i N_i\,\psi^{(j)}\Bigr) ,

.. (vv-status rationale) The governing iteration of the within-group
   splitting, stated so the finalize can be derived from it.  Not a solver
   claim: the map is the DEFINITION of what the SI driver does, and the
   verifiable content — that ONE application at the converged iterate
   reproduces that iterate, and that its angular integral reproduces the
   reported scalar flux — is pinned by
   ``tests/sn/solve/test_eigenvalue_finalize_reconstruction.py``
   (``@pytest.mark.catches("ERR-083")``, seven arms × two orders).
.. vv-status: sn-finalize-map documented

which is the source-iteration bullet above written once.  At a fixed point
:math:`G(\psi^\star) = \psi^\star`, so applying :math:`G` **once** to a
converged iterate returns that iterate — and that is the whole
reconstruction.  Nothing is re-solved, nothing is re-selected: the
finalize reads :math:`M`, :math:`\{N_i\}` and :math:`\psi_{\rm conv}` off
the :class:`~orpheus.sn.solver.InnerSolve` record the last within-group
solve left behind — ``splitting.implicit`` and the ``driven_gains`` the
driver actually applied — and evaluates
:func:`~orpheus.numerics.iteration.fixed_point_step` on them.  That
record carries **both faces of the O-3 query contract** in one object:
``system`` answers with the Problem's operators as POSED, ``splitting``
with the Strategy's operators as USED
(:ref:`sn-splitting-is-a-strategy-value`).

Two things the identity :math:`G(\psi^\star) = \psi^\star` buys are worth
spelling out, because they are what makes ONE step enough rather than a
convenience:

* :math:`M^{-1}` need not be the *representative* the iteration ran.  On
  the windowed arm the driver's step is :math:`P\,M^{-1}` (project to
  moments after the sweep); the record keeps the **un-wrapped** forward
  :math:`M`, so one step through :math:`M^{-1}` alone comes back
  per-ordinate.  A moment iterate is un-windowed by the reconstruction
  *for free*, because the fixed-point identity does not care which
  right-inverse of :math:`M` you use.
* The right-hand side may legitimately differ from the one the iteration
  last saw — which is exactly what the finalize exploits, below.

**The source is the CONVERGED fission source, and that choice is the exit
balance's question.**  The last inner solve ran with
:math:`q = F\phi_{N-1}/k_{N-1}` — the *penultimate* outer's fission source,
because :func:`~orpheus.numerics.eigenvalue.power_iteration` builds
:math:`q` from the iterate it *enters* the outer with and only then
renormalises :math:`\phi` and updates :math:`k`.  The
finalize re-poses with :math:`q_F(\phi_{\rm conv}, k_{\rm conv})`, the
fission source built from the values the caller is actually handed.  The
reason is not tidiness: a caller who checks the returned object is asking
*does this* :math:`\psi` *solve the equation with this* :math:`k` *and this*
:math:`\phi`?, and only the converged source makes the answer yes.  The
shipped exit-balance diagnostic asks exactly that question
(:ref:`sn-exit-balance-projection`), and it would be measuring a
discrepancy of the finalize's own making if the finalize had used the
last inner's source.  At convergence the two sources agree to the outer
tolerance, so the choice costs nothing; at a truncated exit it is the
difference between a diagnostic that tracks the truncation and one that
does not.

**Everything else in the composite arrives as a gain, not as hand-staged
data.**  The lagged couplings are the *value's* own — the
:attr:`~orpheus.sn.splitting.Splitting.explicit` pieces of the
:class:`~orpheus.sn.splitting.Splitting` the inner ran:
:math:`(S, N_{2n}, B_a)` on the Jacobi arm, :math:`(S, N_{2n}, B_{\rm upper})`
when the inner ran under the boundary Gauss-Seidel schedule, and the one
coupled gain grid on a carrying (curvilinear) mesh.  Three consequences:

* **The reflective boundary is** :math:`B\,\psi_{\rm conv}`, delivered
  through ``rhs.boundary`` exactly as in every inner iterate.  The
  finalize does **not** reflect the converged trace into its own inflow
  slots by hand; the inflow is a solved unknown carried on
  ``ψ.boundary`` (Wave O #208 O.4a.2), and the reconstruction re-derives
  it from the same operator the iteration used.
* **The** :math:`\ell \ge 1` **emission is present because it is the
  gain's**, not because the finalize remembered it.  This is the whole
  content of the #448 repair — see the retirement note below.
* **The schedule the inner chose is the schedule the finalize inherits.**
  A boundary-G-S inner reconstructs through :math:`(L+C-B_{\rm lower})^{-1}`
  with :math:`B_{\rm upper}` lagged, which is a *different* splitting of
  the same :math:`A`; the converged answer must not depend on it, and that
  is a pinned row
  (``TestTheGaussSeidelArmPosesItsOwnSplitting`` →
  ``test_the_schedule_does_not_move_the_converged_answer``).

**On a carrying mesh the pair is reconstructed as a pair.**  The
eigenvalue right-hand side is built once, by
:func:`~orpheus.sn.solver._eigenvalue_driver_source` (``[M]`` **3** call
sites — the SI inner, the Krylov inner, and this finalize).  Its System-A
member is the fission source lifted per-ordinate through
:meth:`AngularSourceSink.from_isotropic
<orpheus.transport.source_sinks.AngularSourceSink.from_isotropic>` with a
**zero** external boundary; its System-B member is the :math:`\ell = 0`
fold of the **fission source alone** as the :math:`\psi_{1/2}` march's
entry (#282 route (a)).  The seed folds fission alone because the coupled
gain grid carries the rest — the ``Emission`` and :math:`B_b` blocks — so
folding the total source there would double-count it.  That is not a
finalize-specific rule: it is what the two inner drivers pass, and the
finalize passes it because it calls the same constructor.

.. note::

   ⛔ ``max_outer = 0`` is a legal call that runs the power loop zero
   times, so there is no within-group solve to reconstruct from and
   :func:`~orpheus.sn.solver.solve_sn` **raises**.  This is a live refusal
   of a reachable state, not decoration: the finalize is one step of the
   iteration, and there is no cold-solve fallback by design — a fallback
   would be a second reconstruction path, which is the twin this section
   exists to have retired.

.. warning::

   ⛔ **Until 2026-09-06 this block built its own source, and it was P0
   only** (:doc:`ERR-083 </theory/verification/error_catalog>`).  It
   assembled :math:`F\phi/k + \Sigma_{s,0}^{\mathsf T}\phi +
   \nu_{2n}\Sigma_{2,0}^{\mathsf T}\phi` through three now-retired
   solver-side delegators and lifted it **isotropically**, so at every
   ``scattering_order ≥ 1`` the :math:`\ell \ge 1` half of *both*
   collision channels was absent from the reconstruction's right-hand
   side while the loss arm the iterate converged against carried it.  The
   returned :math:`\psi` therefore solved a different equation from the one
   the solve converged, and its own angular moment did not reproduce the
   :math:`\phi` shipped beside it: ``[M]`` on the 421-group Be-reflected
   slab at :math:`L = 2`, the returned flux missed the converged iterate by
   **8.776e-02** and missed its own reported :math:`\phi` by **3.405e-02**;
   the one-step reconstruction reads **1.236e-10** and **3.170e-10** on the
   same solve.  :math:`k` and :math:`\phi` were never affected — they are
   the power iteration's — which is precisely why every eigenvalue-value
   gate in the tree was structurally blind to it.

   Two things the repair also removed, both worth recording because their
   absence looks like a regression until you check:

   * the **hand reflect** of the converged trace into its own inflow slots
     (``ψ.inflow ← B·ψ.outflow`` before the sweep).  ``[M]`` on a converged
     exit it was INERT — skipping it moved the answer by 2.0e-13 / 2.3e-15
     and bit-identically on a vacuum arm — because the converged inflow
     already equals :math:`B\,\psi_{\rm outflow}`.  No value gate in this
     tree could witness its removal, so the honest artefact is the
     measurement plus a *wrong*-:math:`B` mutation arm, not a gate.  The
     whole-trace verb itself survives as the sweep-tier gates' inter-sweep
     helper (``tests/sn/_test_helpers.py::reflect_outflow_into_inflow``);
     it has **no production caller** any more.
   * the ``AngularBoundarySourceSink.prescribed_inflow`` cast the finalize
     used to perform on that reflected trace (the ERR-071 role conversion).
     The finalize passes no trace at all now — its external boundary source
     is zero and :math:`B` is a gain — so that call site is moot.  The
     factory and the rest of the ERR-071 fix are untouched.

.. _sn-convergence-contract:

The convergence contract — a best-effort answer says so
--------------------------------------------------------

Every entry above can stop for two structurally different reasons: it
**converged**, or its **budget ran out** and the returned field is a
best-effort iterate, mid-descent.  Both come back as the same type from
the same call, so the distinction has to be carried explicitly or it is
lost.

**Where the fact comes from.**  The loop that stops knows why it stopped,
and since 2026-08-08 it says so rather than discarding it.  Since #340 it
says so by **measuring**, never by asserting: each level reports the
quantities it stops on — magnitude and tolerance together — as a
:class:`~orpheus.numerics.convergence.StoppingCriterion`, and every
convergence verdict in the stack is *derived* from those trajectories.

* the **inner** (fixed-source) fact is the driver's own
  :class:`~orpheus.numerics.convergence.IterationRecord`, returned by
  :class:`~orpheus.numerics.iteration.SourceIteration` and
  :class:`~orpheus.numerics.iteration.KrylovAcceleration` alongside the
  iterate;
* the **outer** (eigenvalue) fact is
  :attr:`~orpheus.numerics.eigenvalue.PowerIterationOutcome.converged`, a
  property over the outer record that
  :func:`~orpheus.numerics.eigenvalue.power_iteration` assembles from the
  readings :meth:`SNSolver.measure_stopping_criteria
  <orpheus.sn.solver.SNSolver.measure_stopping_criteria>` returns each
  outer — ``dk`` against ``keff_tol``, ``dphi`` against ``flux_tol``;
* the outer record carries the inner records as **children**, so the two
  facts compose rather than collapsing.

.. important::

   The two questions are genuinely different and the difference is the
   point.  ``record.converged`` asks whether THIS level met its own
   criteria; :attr:`~orpheus.numerics.convergence.IterationRecord.fully_converged`
   asks whether it and every level beneath it did.  An
   increment-only outer stop cannot see an upstream throttle — a truncated
   inner suppresses the very increments the outer reads, so the outer
   stalls and calls the stall convergence.  `[M]` on a 20-cell 2-group
   :math:`S_8` slab at ``max_inner=1`` the outer reports
   ``converged=True`` with :math:`\keff` wrong by **11×** its own
   ``keff_tol``; ``fully_converged`` is ``False`` and
   :attr:`~orpheus.numerics.convergence.IterationRecord.first_failure`
   names the starved inner.  **A value gate asserting physics must read
   the fold, not the level.**

⛔ Until 2026-08-09 the inner fact was "one shared predicate over the
residual history" (``orpheus.sn.solver._claims_convergence``).  That
predicate is retired: it read an EMPTY history as *not converged*, so a
Krylov solve that returned on its initial guess was indistinguishable from
a truncation, and reusing it as an audit instrument produced `[M]` 44 of 90
phantom truncations.  A record separates the two —
:attr:`~orpheus.numerics.convergence.IterationRecord.iterated` is the
discriminator.

These surface directly on the Solution:
:attr:`~orpheus.sn.solution.SolutionBase.record` **is** the tree, and
:meth:`SolutionBase.converged() <orpheus.sn.solution.SolutionBase.converged>`
reads its
:attr:`~orpheus.numerics.convergence.IterationRecord.fully_converged`
fold.  Nothing sits between the two.  Every diagnostic is asked of the
object that owns it — the **record** for the path, the **outcome** for the
answer, the **certificate** for what the exit measured about the returned
state — and the next section says which reading moved where, and which
one deliberately did not move at all.

⛔ Until 2026-08-09 ``converged`` was a *field*, and this paragraph argued
its honesty from the fact that it was **required**: no default, so a
producer could not claim convergence by omission.  That was the right fix
for the wrong layer — a required field still has to be WRITTEN, by hand,
at every producer, and #342 was five such writes with one of them a literal
``True``.  Deriving it removes the question of who writes it: there is no
argument to pass, so there is nothing to get wrong.

.. tip::

   Read :attr:`Solution.record
   <orpheus.sn.solution.SolutionBase.record>` for anything the flat
   readings drop — which is most of it.  ``record.report()`` prints the
   whole tree, level by level, with each criterion's last value, the
   tolerance it was judged against, the observed rate, and the budget that
   rate projects; it is written to be pasted into a bug report unedited.

.. warning::

   **A value gate that does not assert convergence is asserting an
   arbitrary iterate.**  This is not hypothetical.  `[M]`
   ``test_d3_pure_absorber_per_ordinate_psi_exact`` asserted a closed-form
   identity to ``rtol=1e-10`` on an all-reflective 3-D box that needs
   **1631** sweeps, against the **then-default** ``max_inner`` of **1000**
   (hardcoded; the default has been derived from the tolerance since #340
   N3 landed 2026-08-09).  It read
   the 999th iterate, never read the flag the solver had honestly set to
   ``False``, and passed for months because the truncated error happened to
   land inside the tolerance — until a *correct* quadrature change (#337)
   moved it out.  The one-line defence is to assert
   ``sol.converged()`` **before** reading any value — it reads the
   TREE-wide predicate, because on an eigenvalue solve a LEVEL-wide
   ``converged`` reads ``True`` while an inner starves (see *Loudness*
   below).

   The diagnostic tell is worth memorising, because it points the wrong
   way: the error was **bit-identical at every** ``inner_tol`` **from 1e-9
   to 1e-15**, which reads as a discretization floor.  It is the opposite —
   the running residual never fell below even the loosest tolerance, so
   every run hit the same cap and returned the same bytes.
   Tolerance-insensitivity means the tolerance never *bound*; read the
   iteration count against the budget before concluding anything about the
   discretization.

   ⚠ Read it against
   :attr:`~orpheus.numerics.convergence.IterationBudget.in_iterations`, not
   against the raw knob you passed. The two coincide for source iteration,
   the power outer, CP and MoC, and they do **not** for Krylov: ``max_inner``
   there is scipy's restart-CYCLE cap, while the recorded trajectory counts
   inner Arnoldi steps, and one cycle buys ``restart`` of them (which
   :func:`~orpheus.sn.solver` sizes to the full ``n_dof``). Comparing the raw
   pair is ERR-079 (#349) — it read a healthy converged solve as
   having run out. :meth:`~orpheus.numerics.convergence.IterationRecord.report`
   already prints the honest ceiling, so prefer it to hand arithmetic.

**Loudness.**  A truncated exit emits
:class:`~orpheus.numerics.convergence.ConvergenceWarning`, naming **the level
that failed** and the budget *that level* ran out of, the tolerance *its*
binding criterion missed, and how far its last iterate was — the distance
between "one more sweep" and "diverging" wants opposite responses.  The level
matters because a solve is a tree: on an eigenvalue run it is usually the
*inner* that starved, and until #340 N6 the message quoted the entry's own
``max_outer`` instead — a knob that cannot help, since the outer's stop test
is entirely increments and the starved inner is what suppresses them.

The guard is
:attr:`~orpheus.numerics.convergence.IterationRecord.fully_converged` —
**every** level, not the top one — so a converged outer standing on a starved
inner is audible.  That case is the whole of the #340 headline defect, and it is not
exotic: `[M]` 2026-08-10, **20 tests in the shipped suite** sat in exactly
that state, at observed rates :math:`\rho` between 0.889 and 0.993, every one
of them silent.

.. note::

   **Ask** :meth:`sol.converged() <orpheus.sn.solution.SolutionBase.converged>`
   — which is
   :attr:`~orpheus.numerics.convergence.IterationRecord.fully_converged`, the
   FOLD — **not** the record's per-level ``converged``, before asserting
   physics against a result.  The two differ precisely on the starved-inner
   solve: the level-wide reading is ``True`` there, because the outer really
   did meet its own criteria — it just met them on increments an upstream
   throttle had suppressed.

⛔ Until 2026-08-10 the guard was ``converged``, the TOP level only, and the
widening was scheduled to ride an *outer residual certificate* that would
separate a truncation which corrupted the answer from one that did not.  `[M]`
that certificate was **refuted by measurement**: the benign and corrupting
populations overlap **634×** and it misses 15 of 16 corrupting cases.  The
ruling that followed was to widen the guard *unconditionally* — a truncation
the caller has not declared is worth saying out loud whether or not we can yet
say what it cost.  The 20 silent tests were adjudicated in the same change —
10 declare the truncation as their fixture and suppress this one category
in-test, 10 are audible on purpose and tracked with measured budgets in
`#352 <https://github.com/deOliveira-R/ORPHEUS/issues/352>`_.

.. note::

   **This machinery is no longer SN's** (#340 N4.7, 2026-08-11).  The emitter
   left ``sn/solver.py`` for
   :func:`~orpheus.numerics.convergence.warn_if_unconverged`, and CP, MoC and
   1-D diffusion now call it from their own public entries.  Nothing about
   SN's behaviour changed — `[M]` the emitted message is character-identical
   across all four advice arms plus the nested and balance-defect cases,
   verified against the pre-move function lifted out of git — but two facts
   about the SHAPE of the diagnostic are worth carrying:

   * The helper was already ~90 % family-agnostic. Every fact it reads off the
     failing level is a generic
     :class:`~orpheus.numerics.convergence.IterationRecord` member; only
     ``balance_defect`` was SN's, and it is now an optional keyword —
     typed :class:`~orpheus.numerics.outcome.Evidence` since 2026-09-17 —
     that the other three leave at its ``NotApplicable`` default.  Only a
     :class:`~orpheus.numerics.outcome.Measured` value renders; every
     other member renders an *absent* clause, never the word
     "unavailable", because an empty clause cannot be misread as a
     measurement.
   * ⛔ Its closing advice used to name the literal string
     ``solution.history.fully_converged`` (``[M]``
     ``orpheus/sn/solver.py:594`` at ``28435e11``; the flat ``history`` view
     that string named is itself retired —
     :ref:`sn-the-record-answers-for-its-own-level`).  That is a guess at the
     CALLER's local variable name — a fact no library can know — and it was
     outright wrong for the three families whose entries return a
     ``*Result``.  It now names the attribute and its type.  A per-entry spelling passed in as an
     argument was considered and **rejected**: it would re-commit the exact
     defect N6a retired, a fact asserted by the call site and free to drift
     from the object it describes.

.. _sn-the-record-answers-for-its-own-level:

What the record answers, and what deliberately did NOT move onto it
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

⛔ **Until step 3's unit U6 (2026-09-17) a flat ``IterationHistory``
dataclass sat between the Solution and its record**, and the two steps of
its life are both worth keeping.  It was born as six hand-written FIELDS
(Issue #197 PR-TYPED-5), which is how ``n_inner`` came to mean
``len(residuals)`` on one path and ``len(residuals) + 1`` on the other —
undocumented and exactly backwards from the truth — and how ``converged``
came to be written five times by hand with one of the five a literal
``True`` (#342).  A projection maintained by hand at :math:`N` sites is
:math:`N` chances to disagree.  From 2026-08-09 (#340 N2b-ii) it became a
**view**: every scalar DERIVED from the record, so the flat surface could
no longer drift from the tree it summarised.  That was the right repair,
and it left a second defect standing — a *second surface* naming the same
facts, which a reader had to choose between, and onto which no new
question could be added without widening the very flattening the record
exists to undo.  Its own docstring said so — *"do not grow this surface:
add the question to the record, where the tree can answer it"* — and U6 is
that instruction carried out.

.. list-table:: the retired view's readings, and what answers now
   :header-rows: 1
   :widths: 30 70

   * - retired reading
     - what answers it now
   * - ``converged`` / ``fully_converged``
     - :attr:`record.converged
       <orpheus.numerics.convergence.IterationRecord.converged>` /
       :attr:`record.fully_converged
       <orpheus.numerics.convergence.IterationRecord.fully_converged>` —
       the same two properties the view delegated to, one hop shorter
   * - ``keff_history``
     - :attr:`outcome.trajectory
       <orpheus.numerics.outcome.EigenOutcome.trajectory>`.  It was never
       a convergence quantity: what the outer *stops* on is the
       per-iteration INCREMENT ``dk``, which lives in the record; the
       :math:`k` sequence is a **physics output** and belongs on the
       answer, co-indexed with :math:`\lambda`
       (``trajectory[-1] == lam``)
   * - ``dominance_ratio()`` / ``latest_keff()``
     - :meth:`outcome.dominance_ratio()
       <orpheus.numerics.outcome.EigenOutcome.dominance_ratio>` and
       ``outcome.lam`` — readings of that trajectory, so they live beside
       it.  On a ``Solution[SourceOutcome]`` neither **exists**, which is
       the type saying what a ``None`` used to say badly
   * - ``balance_defect`` / ``gauge_correction``
     - :attr:`certificate.balance
       <orpheus.numerics.outcome.ExitCertificate.balance>` /
       :attr:`certificate.gauge
       <orpheus.numerics.outcome.ExitCertificate.gauge>`, as typed
       :class:`~orpheus.numerics.outcome.Evidence`.  These were the two
       ``float | None`` members carrying five and three documented
       meanings; the sum names the reason instead of erasing it
       (:ref:`the-solution-outcome`)
   * - ``flux_residuals`` / ``latest_residual()``
     - :attr:`record.trajectory
       <orpheus.numerics.convergence.IterationRecord.trajectory>` — NEW
       at U6; see below
   * - ``total_inner_iterations``
     - :attr:`record.leaf_iterations
       <orpheus.numerics.convergence.IterationRecord.leaf_iterations>` —
       NEW at U6; see below
   * - ``n_inner`` / ``n_outer``
     - :attr:`record.n_iterations
       <orpheus.numerics.convergence.IterationRecord.n_iterations>` of the
       solve's own top record.  **Deliberately not two names**; the
       argument is the subsection below
   * - ``_is_outer``
     - nothing.  It existed only to let one name serve two kinds, and the
       kind is the Solution's TYPE now.  (It read ``bool(record.children)``
       — the tree's own structure, never the level's ``label``, because
       reading a control decision off a string chosen for humans is
       stringly-typed dispatch.)

**Two readings the record gained**, because they are questions about a
tree and the tree is the thing that can answer them:

* :attr:`~orpheus.numerics.convergence.IterationRecord.trajectory` — the
  **binding criterion's** per-iteration trajectory, ``()`` when nothing
  bound.  For a within-group solve that is the relative flux increment,
  which is exactly what ``flux_residuals`` meant; ``()`` is what a
  consumer asking *"did this level iterate on a residual?"* needs to read
  when it did not, and the DSA rate diagnostics branch on precisely that,
  so their ``if not …trajectory:`` guards keep their meaning verbatim.
* :attr:`~orpheus.numerics.convergence.IterationRecord.leaf_iterations` —
  the iterations run by the **leaves** of this tree, which is the *inner
  work*: a leaf reports its own count, an outer reports the sum over every
  leaf beneath it.  It is the measurand of the SI spectral-rate and
  Gauss-Seidel-recovery diagnostics — it, not the outer count, is what a
  preconditioner is meant to lower — and it is never ``None``: a level
  that never iterated reports ``0``.

.. note::

   ``leaf_iterations`` is a **generalisation** of the reading it replaces,
   not a re-spelling.  The retired ``total_inner_iterations`` summed
   ``n_iterations`` over the record's DIRECT children; ``leaf_iterations``
   recurses to the leaves.  On a tree deeper than two levels those are
   different numbers, and the retired form reports the intermediate levels'
   own counts while never reaching a leaf: ``[M]`` 2026-09-17, on a
   hand-built ``outer(3) → inner(5) → {gmres(11), gmres(13)}`` record the
   retired sum reads **5** and ``leaf_iterations`` reads **24**.

   On the shipped S\ :sub:`N` path the two agree, because the eigenvalue
   record is exactly **two** levels deep — ``[M]`` 2026-09-17, ``A-2g``
   10-cell slab, ``gauss_legendre(8)``, at ``solve_sn``'s default
   tolerances: ``outer(power-iteration)`` with three LEAF children, retired
   sum **394** and ``leaf_iterations`` **394** under source iteration (both
   ``gauss_seidel`` and ``jacobi``), and **256** / **256** under Krylov
   (``inner(gmres)``).  So U6 is bit-identical on every reading this page
   reports, and the generalisation is insurance against the day a third
   level appears.  The deep case is pinned on the record itself, where it
   belongs:
   ``tests/numerics/test_iteration_record.py::TestTheTreeReadingsAreTheRecords``.

The ``None``-by-SHAPE trio did not move, and that is the point
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Three of the view's readings were typed ``int | None`` — ``n_inner``,
``n_outer`` and ``total_inner_iterations`` — and where the ``None`` was
reachable it did not mean *"not measured"*.  It meant **"you asked the
wrong kind"**: ``n_inner`` read ``None`` whenever the top level had
children (an eigenvalue solve), ``n_outer`` whenever it did not (a
fixed-source one).  ``flux_residuals`` carried the same convention in the
other direction, returning ``()`` on an OUTER level so that one name could
be taken on either kind.

⚠ The third member's ``| None`` was never reachable at all — both branches
of ``total_inner_iterations`` returned an ``int``, as its own docstring
said (*"populated by both paths"*).  ``[M]`` 2026-09-17, the retired
property lifted verbatim out of ``git`` and run against a live record
(``A-2g`` 10-cell slab, GL-8, source iteration): on the outer it reads
``394`` where ``n_inner`` reads ``None``, and on that outer's leaf inner it
reads ``190`` where ``n_outer`` reads ``None``.  So the trio was two
kind-keyed ``None``\ s plus one piece of annotation debt that the flat
surface's Optional-by-shape habit had spread onto a total that was always
a number — which is the shape of the defect, not a detail of it: a
convention adopted for two members is copied onto the third for symmetry,
and every consumer then writes a guard for a state that cannot occur.

Those ``None``\ s existed **so that ONE name could be read on either
problem kind**.  Since step 3 U2 the kind is the Solution's TYPE
(``Solution[EigenOutcome]`` versus ``Solution[SourceOutcome]``), so the
question the ``None`` answered can no longer be asked — and re-minting the
trio on :class:`~orpheus.numerics.convergence.IterationRecord` under the
old names would have carried a *kind* discriminator onto an object that
has no kind.  A record is a level in a tree; ``n_inner`` and ``n_outer``
are not two quantities it holds, they are one quantity — the iterations
THIS level ran — wearing two names chosen by what the level happened to be
nested inside.

The census says the distinction was never load-bearing at a single call
site.  ``[M]`` 2026-09-17, by the U6 diff (every read the carve migrated,
alias receivers included; predicate = an attribute read of the named
member on a ``.history`` receiver or on a local bound to one):

.. list-table:: every migrated reader of the trio and of ``flux_residuals``
   :header-rows: 1
   :widths: 26 10 64

   * - reading
     - reads
     - where, and what kind of solve each one holds
   * - ``n_inner``
     - 23
     - ``test_dsa_acceleration`` 2, ``test_dsa_rate`` 2,
       ``test_si_convergence_rate`` 6, ``diag_d3_absorber_01`` 1,
       ``diag_d3_absorber_02`` 9, one ``_solve_fixed_source_si`` docstring,
       plus 2 rows of the view's own retired tests — **every consumer holds
       a ``solve_sn_fixed_source`` result**
   * - ``n_outer``
     - 10
     - ``diag_vacuum_bc_eigenvalue_divergence`` 3,
       ``test_sn_adjoint_entries`` 2, ``test_si_convergence_rate`` 2, plus
       3 rows of the view's own retired tests — **every consumer holds an
       eigenvalue (or adjoint eigenvalue) result**
   * - ``total_inner_iterations``
     - 4
     - ``test_si_convergence_rate`` 1 (the eigenvalue-path measurand) plus
       3 of the view's own retired rows
   * - ``flux_residuals``
     - 11
     - ``test_dsa_acceleration`` 4, ``test_dsa_rate`` 2,
       ``diag_d3_absorber_01`` 1, ``test_convergence_contract`` 1, plus 3
       of the view's own retired rows — all fixed-source, all reading the
       within-group increment

So both names were always ``record.n_iterations`` of the solve's own top
record, read by a consumer that already knew which kind it was holding —
and the ``None`` branch those call sites carried (``if history is None or
history.n_inner is None: pytest.fail(…)``) was guarding against a state
its own fixture made unreachable.  ``[M]`` the carve removed **four** such
consumer guards — two on ``n_inner`` (the DSA rate and acceleration
diagnostics), one on ``n_outer`` (the adjoint-entries gate) and one on
``total_inner_iterations`` (the SI-rate gate, whose failure message named
a state that could not occur) — plus the **two** assertions in the view's
own tests that pinned the convention itself (``h.n_inner is None`` on an
outer, ``h.n_outer is None`` on a leaf).  Those last two are the tell: the
only rows that ever exercised the ``None`` were the ones written to
document it.

``record.trajectory`` widens the ``flux_residuals`` reading in exactly one
way, and it is the honest one: **an outer level now answers with its own
binding criterion's trajectory** — the eigenvalue increments — where the
view returned ``()`` so the name could be shared.  Nothing in the shipped
tree read ``flux_residuals`` on an eigenvalue Solution — the table above
is the denominator, and its only eigen-path rows are the view's own tests
— so no consumer changed meaning; a future one asking an
outer for its trajectory gets the answer for its own level instead of an
empty tuple it must know to distrust.


.. _sn-exit-balance-projection:

The balance projection
----------------------

The refutation left a real question standing — *how much did this truncation
cost?* — and the answer the warning carries is the **per-group neutron-balance
defect of the returned iterate**, reported on the Solution's exit certificate
as :attr:`ExitCertificate.balance
<orpheus.numerics.outcome.ExitCertificate.balance>`:

.. math::
   :label: sn-exit-balance-defect

   R_g \;=\; \int_V \int_{4\pi} \bigl(A\psi - q\bigr)_g \, d\Omega \, dV,
   \qquad
   \text{reported as } \;\frac{\lVert R_g \rVert}{\lVert R_g(q) \rVert}

**Why the projection and not the residual norm** — and the reason is
structural, not statistical.  Up to **99.995 %** of :math:`\lVert r \rVert` is
reflective-trace rows, and a reflective inflow-trace defect in a zero-leakage
system carries **no net current**, so a balance-based :math:`\keff` is blind
to it *by conservation*; `[M]` the transfer gain
:math:`\lvert\Delta k\rvert / \text{defect}` spans **1.16 × 10⁵**.  Integrating
over angle and volume annihilates exactly those rows, because that is the
functional :math:`\keff` itself reads.  `[M]` the overlap falls **634× →
4.64×**.

.. warning::

   **4.64× is still an overlap.  This is a diagnostic magnitude, never a
   threshold.**  Do not branch on it, and do not assert on it in a test — a
   gate that did would be the refuted certificate wearing a different name.
   It is reported so a reader who has been told their solve truncated can
   weigh how much that is likely to have cost.

   ⛔ Nor can it be sharpened with a cheap adjoint: `[M]` a spatially-flat
   0-D weight makes it **worse**, 4.64× → **128.95×**, because a signed
   projection against a wrong weight manufactures near-cancellations, i.e.
   false negatives.  The weighting channel already exists
   (:meth:`IntegratedReactionRate.evaluate` takes ``adjoint=``); what a real
   gate needs is the adjoint *solve*
   (`#350 <https://github.com/deOliveira-R/ORPHEUS/issues/350>`_).

It is computed **only on the exit that warns** — the exact complement of the
within-group certificate, which fires when the solve *claimed* convergence and
raises.  One equation, two verbs, complementary guards, so no solve pays for
both forward applies and the converged path costs what it always did.  `[M]`
one residual evaluation is ≈ 3 inner iterations, i.e. **0.72 %** of a
400-iteration truncated solve.

.. warning::

   ⛔ **The complement of a guard is not the same as coverage, and this
   pair left the exit uncovered until 2026-09-06** (#448 /
   :doc:`ERR-083 </theory/verification/error_catalog>`).  The two verbs
   above are complementary *in when they fire* — certificate on the
   converged exit, defect on the truncated one — and they are NOT
   complementary in *what they read*.  The certificate reads the within-group
   **iterate**, inside the inner solves; this projection reads the returned
   flux, but only when the solve did not converge.  So the object a caller
   receives from a **converged** eigenvalue solve was evaluated by neither,
   and the hand-built P0-only reconstruction that produced it drifted
   undetected for the whole life of the anisotropic solve.

   Both halves are now repaired.  The finalize is one step of the
   iteration's own map (:ref:`sn-finalize-one-step`), so the returned flux
   solves the equation the reported :math:`(k, \phi)` pose; and the
   projection's budget response is the property a regression gate asserts
   (``tests/sn/solve/test_eigenvalue_finalize_reconstruction.py`` →
   ``TestTheShippedDiagnostic``).  ``[M]`` 2026-09-06, 2-group A|B|A slab, ``keff_tol =
   flux_tol = 1e-12``, ``inner_tol = 1e-11``, ``max_outer`` 3 → 12: the
   defect falls by :math:`1.43\times10^{7}` at :math:`L = 0` and
   :math:`3.46\times10^{7}` at :math:`L = 1`.  Before the fix the
   :math:`L = 1` column fell by **1.0002 ×** — pinned at a floor set by the
   reconstruction rather than by the truncation, i.e. a diagnostic that
   could not move is a diagnostic that carried no information
   (``vv-principles`` #19).  The full before/after table is in the ERR-083
   entry.

   ⚠ Reading note for anyone comparing the two trees: the absolute defect
   at a given budget is **not** comparable across the fix, because the two
   finalizes construct the returned flux differently — the pre-fix one
   re-solved a source built from :math:`\phi`, the post-fix one steps the
   map from the iterate, and at a truncated exit those differ by the
   truncation itself (compounded by the ERR-052 between-solve
   renormalisation).  What is comparable, and what the diagnostic
   advertises, is the RATIO down the budget.

**When no number is reported, the certificate says WHY** — and that is
the second thing step 3 changed here.  Until 2026-09-17 the answer was a
``None`` carrying five documented meanings, which is stringly-typed
dispatch wearing an absence: the consumer had to re-derive which one
applied from context the type did not carry.  The member is now a typed
:class:`~orpheus.numerics.outcome.Evidence` sum
(:ref:`the-solution-outcome`), and the SN evaluator produces every branch
of it:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - what comes back
     - when
   * - :class:`~orpheus.numerics.outcome.Measured`
     - the exit that warns — a truncated solve, the case this whole
       section is about
   * - :class:`~orpheus.numerics.outcome.Certified`
     - the tree **fully converged**: the within-group exit certificate
       already ASSERTED :math:`\lVert A\psi - q\rVert/\lVert q\rVert`
       within its safety factor (raising otherwise), so no forward apply
       is spent and the *bound* is reported instead of a number
   * - :class:`~orpheus.numerics.outcome.NotApplicable`
     - the source integrates to zero per group — the ratio is
       **undefined**, which is a different statement from *unmeasured*
   * - ``NotYet(310)``
     - **moment-tailed (LD) schemes** on both fixed-source arms: the
       residual mint does not admit the trailing :math:`2^d`
       spatial-moment axis, so no residual exists to project.  The same
       un-built widening the within-group certificate exempts, reached
       through one shared predicate rather than two copies of the test
   * - ``NotYet(353)``
     - **the daggered eigenvalue entry** ``solve_sn_adjoint``: N5's
       population is forward-only, so there is no reference against which
       a plausible number could be checked.  Deferred deliberately rather
       than guessed — assembling one from plausibility is the ERR-032
       class

⭐ **One of the two old** ``None`` **entries came back.**  The
**carrying** (curvilinear) eigen arm used to report nothing, and the
reason was a real refusal rather than an oversight: what
:func:`~orpheus.sn.solver.solve_sn` assembled at its exit was a *bare*
System-A residual against a System-A fission rhs, which on a carrying
mesh silently omits :math:`r_B` — the ``vv-principles`` Mode-12 blindness
the split-residual mint exists to prevent — so "no number" was more
honest than a residual missing a block.  Since the certificate reads the
**outcome's own** residual and rhs, and the outcome's question is the
coupled pencil, the missing piece assembles itself: the rhs is
:math:`\mu(\lambda)\,M\psi` on the coupled space, and
`#354 <https://github.com/deOliveira-R/ORPHEUS/issues/354>`_ — whose gap
*was* the un-assembled coupled rhs — is measurable.

⚠ The mirror case is a **repair, not a re-baseline**, and it moves a
published number.  :func:`~orpheus.sn.solver.solve_sn_multiplying_source`
solves :math:`(A - F)\psi = q` and its defect was being computed against
:math:`A` alone — the residual of a pure-transport problem that entry
does not pose.  ``[M]`` 2026-09-17 on the witness's truncated subcritical
slab (:math:`L = 4`, hub :math:`k = 0.907457573`, ``inner_tol = 1e-12``,
``max_inner = 5``): :math:`0.8294593510371534` against :math:`A`
:math:`\to` :math:`0.8758249879057027` against :math:`A - F`.  The other
four entries pose ``pencil.at(0)``, which **is** the loss operator by
object identity, so their numbers did not move at all.

It is a warning rather than an exception by the ERR-053
precedent (legitimate callers harvest the residual history of a
deliberately-truncated solve), and it escalates to a hard failure with::

    python -O -m pytest -W error::orpheus.numerics.convergence.ConvergenceWarning

.. warning::

   ⛔ **The category must be DOTTED, and this page said otherwise until
   2026-08-09.**  The recipe was published as ``-W
   error::ConvergenceWarning`` here and at four code sites (one of them the
   emitted warning message itself).  That string **does not parse** —
   Python resolves an undotted ``-W`` category against ``builtins``, so
   pytest exits at startup with ``AttributeError: module 'builtins' has no
   attribute 'ConvergenceWarning'`` and collects **zero** tests.  The CI
   contract was imaginary for exactly as long as it was documented.

   It survived because the gate that appeared to prove it,
   ``test_it_is_escalatable_to_an_error``, installs the filter through
   ``warnings.simplefilter`` — a true claim about the *category* that says
   nothing about the *spelling*.  The doc, the runtime message and the test
   all agreed on a string no interpreter accepts.

   The spelling is now derived from the class as
   :data:`~orpheus.numerics.convergence.ESCALATION_FLAG`, and
   ``test_the_published_escalation_flag_actually_parses`` consumes that
   STRING through pytest's own parser.  General rule this earned: **for
   every recipe a doc publishes as a command, one gate must consume the
   string, not the API.**

**Budget sizing.**  ``max_inner`` is **derived from the tolerance**, and it
had to be: `[M]` an all-reflective box needs ~an order of magnitude more
sweeps per added dimension (d=1 **32**, d=2 **258**, d=3 **1631**), and the
cost scales as :math:`\Sigma_t \cdot n_{\rm inner} \approx` constant.  One
vacuum face collapses the d=3 figure to 208 — the expensive corner is
specifically zero-leakage, weakly-absorbing, and 3-D.  **A constant cannot
track a tolerance it does not know about**, so ``max_inner=None`` — the
default at every SN entry — resolves through
:func:`~orpheus.numerics.convergence.resolve_iteration_budget` to
:func:`~orpheus.numerics.convergence.default_iteration_budget`, which
inverts the geometric budget law at a stated served rate.  An explicit
``int`` is still honoured untouched: a caller who knows their spectral
radius, or who is deliberately starving a solve to measure its truncated
exit, is exercising the API correctly.

⛔ *Until #340 N3 landed (2026-08-09) the default WAS a hardcoded constant*
— five of them in SN plus a sixth in ``KEigenvalue`` — and `[M]` **both**
SN families were short at d=3 zero-leakage: 830 sweeps needed at
``inner_tol=1e-8`` against 200 shipped, 1441 at ``1e-12`` against 1000.
The two families also differed by **5×** where the law puts the factor at
:math:`\ln(10^{-12})/\ln(10^{-8}) = 1.5`; `[M]` the shipped ratio is now
``1308 : 1961 = 1.4992``.

Gates: ``tests/sn/solve/test_convergence_contract.py`` — each honesty claim
is a PAIR (a converging configuration and a deliberately-starved one),
because asserting ``converged is True`` on a solve that converges is
satisfied by the very hardcoded ``True`` the contract forbids; it is the
starved leg that has teeth.  `[M]` a 6-mutation battery, positive control
first, reddens as designed — including re-introducing #342 verbatim.

.. _sn-exit-gauge:

The exit gauge — a converged solve can still be one of many
------------------------------------------------------------

The convergence contract answers *"did the iteration finish?"*.  There
is a second, structurally different way for a returned field to be
unsatisfying, and no convergence certificate can see it: the **equation**
may not have a unique answer.  On a closed reflective Cartesian box under
diamond differencing :math:`A = L+C-S-N_{2n}-B` is **exactly
singular**, so
:math:`A\psi = q` has a solution *manifold* and the iteration freezes an
arbitrary member of it.  The derivation, the counting laws, the evidence
and the remedy hierarchy are at :ref:`sn-loss-kernel-gauge`; what belongs
here is the exit behaviour.

Every entry that returns a trace applies the :math:`G`-orthogonal
projection :math:`\psi \mapsto \psi - \Pi\psi`
(:eq:`sn-loss-kernel-gauge-projection`) and records the magnitude it
removed on :attr:`ExitCertificate.gauge
<orpheus.numerics.outcome.ExitCertificate.gauge>`.  The
gauge is the **sibling** of the balance projection with one sharpening
that changes what verification it owes: :eq:`sn-exit-balance-defect`
*reports*, and this one *mutates*.  A forgotten balance-defect site
loses a diagnostic; a forgotten gauge site silently returns a
non-physical answer.  Coverage is therefore gated by an enumeration
**derived from the module** rather than hand-listed
(``tests/sn/solve/test_every_entry_gauges_its_trace.py``).

Three properties make firing it at a converged exit safe, and each is
asserted rather than assumed:

* **Residual-neutral.**  :math:`A(\psi - \Pi\psi) = A\psi`, so **no
  convergence certificate can move**.  ``[M]`` on a deliberately
  truncated SI solve the balance defect reads ``0.3111434602740818``
  before and after, while the correction goes
  :math:`3.59\times10^{-2} \to 4.9\times10^{-17}`.  It is applied
  **after** the defect is measured, so the reported number describes the
  object the caller receives.
* **Bulk-invariant.**  The kernel is pure-trace (``[M]`` bulk share
  :math:`1.1\times10^{-28}`), so :math:`\keff`, the scalar flux and
  every reaction rate are untouched — ``[M]`` :math:`\keff` reproduces
  the analytic :math:`\kinf = 1.875` on every mesh, gauged or not.
* **Not a universal absorber.**  ``jacobi`` already lands on the
  canonical member, and there the gauge must remove *nothing*: ``[M]``
  :math:`\sim10^{-15}`, and the pre-gauge deviation measures
  :math:`1.0000` **out of** span.

The gauge member follows the balance member's discipline exactly, and
since 2026-09-17 the type carries it rather than the prose: a
:class:`~orpheus.numerics.outcome.Measured`
:math:`\sim10^{-15}` is the statement *the freedom is real and the solve
landed on the canonical member anyway*, and the two ways of measuring
nothing are **two different values** —
``NotApplicable("no kernel freedom: …")`` when the configuration has none,
and ``NotApplicable("the closure is unclassifiable, so the trace was NOT
gauged: …")`` when the face-mode damping could not be classified.  ⛔
Until then both were the same ``None``, alongside a measured zero-ish
number, on one ``float | None`` field — three states on a type that can
express two.  The third state is the one a caller must never collapse
into the first: an unclassified closure means the trace was **not
repaired**, not that there was nothing to repair.

.. warning::

   :class:`~orpheus.sn.operators.loss_kernel_gauge.GaugeFreedomWarning`
   is **deliberately not** a
   :class:`~orpheus.numerics.convergence.ConvergenceWarning`, and the
   distinction is not cosmetic.  That family means *"an iterative solve
   exhausted its budget; the answer is best-effort"*.  This is the
   opposite situation: ``[M]`` the configuration where it fires hardest
   reports ``fully_converged = True`` and a ``Certified`` balance member
   (the within-group exit certificate asserted the bound, so no number
   was owed).  The solve is fine; the **equation** is degenerate.  Reusing the
   category would also make every caller who escalates
   :data:`~orpheus.numerics.convergence.ESCALATION_FLAG` start failing
   on an unrelated condition.

   It reports an **action taken**, not a configuration property — which
   is what keeps it off the standard all-reflective :math:`\kinf`
   lattice, where it would otherwise fire on every solve.  It escalates
   with :data:`~orpheus.sn.operators.loss_kernel_gauge.GAUGE_ESCALATION_FLAG`,
   whose value is **derived from the class** rather than retyped: the
   category must be DOTTED, for the reason the balance-projection
   section records above.

   The third state is the one a caller must not collapse: an
   **UNDETERMINED** closure — one whose face-mode damping could not be
   classified — warns loudly and is **not** gauged.  ``[M]``
   ``linear_discontinuous`` at :math:`d=3` is exactly that.

.. _sn-solution-carries-its-posing:

What the Solution carries — the answer fused with its question
---------------------------------------------------------------

Everything above describes things the exit *measures*.  This section is
about the object those measurements travel on, and about a defect that
was structural rather than numerical: until 2026-09-17 a
:class:`~orpheus.sn.solution.Solution` recorded the **answer** and not the
**question**, so two solves over one hub could be indistinguishable by
their data while having solved different equations.

The tier-agnostic theory — why the outcome is FUSED, what a gauge *is*
(a section of a torsor quotient), and why the certificate is a sum rather
than a nullable float — is
:ref:`the-solution-outcome`, and the S\ :sub:`N`-specific realization
table is :ref:`the-outcome-sn-realization`.  What belongs here is the
S\ :sub:`N` shape and the consequences for a caller.

Five members, and the KIND is one of their types
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:class:`~orpheus.sn.solution.SolutionBase` is generic in the outcome —
``SolutionBase[O]``, ``O`` **constrained** to
:class:`~orpheus.numerics.outcome.EigenOutcome` /
:class:`~orpheus.numerics.outcome.SourceOutcome` — and carries exactly:

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - member
     - what it is
   * - ``mesh``
     - the **Problem** (the hub), the base point everything else is
       relative to
   * - ``outcome``
     - the kind-typed **answer**, fused with the question it answered,
       the returned STATE and the gauge that picked the representative.
       *The kind IS this member's type*
   * - ``strategy``
     - the :class:`~orpheus.sn.splitting.Splitting` VALUE the solve drove
       — the labelled piece set and the schedule.  Budgets and tolerances
       ride the record, per level, because they are per level
   * - ``certificate``
     - what the exit MEASURED about the returned state, member by member,
       as typed :class:`~orpheus.numerics.outcome.Evidence`
   * - ``record``
     - the Strategy's **path** — the
       :class:`~orpheus.numerics.convergence.IterationRecord` tree

The ROLE stays what #276 A5 made it — a class
(:class:`~orpheus.sn.solution.Solution` /
:class:`~orpheus.sn.solution.AdjointSolution`) — because the verb set
varies by role and **not** by kind.  So the family is two leaves and one
parameter rather than four classes, and the two axes remain what they
always were: role = type, kind = a **type parameter** (it was a
*property* until this step).

⛔ **Until step 3 the kind was read off the ANSWER**: ``keff is not
None``, through ``is_eigenvalue()`` / ``is_fixed_source()``.  Three
things followed, each of which the shape retires rather than documents:

* a **multiplying-source** Solution was indistinguishable from a
  pure-transport one by its data, because both carry ``keff = None`` —
  the two entries solve different equations over the same hub;
* the **eigen gauge was recorded nowhere**, so a consumer holding a
  converged flux could not determine the scale it was on, and two
  perfectly truthful "normalized to unit production rate" states could
  differ by the :math:`(n,2n)` channel;
* ``compare`` branched on ``keff is not None``, so a cross-kind pair
  silently skipped the eigenvalue channel instead of refusing.  It now
  refuses on **three** axes — role, kind and phase space — and its
  ``keff_abs`` channel is ``Evidence`` rather than a nullable float.

``keff``, ``is_eigenvalue()``, ``is_fixed_source()``, ``keff_history``,
``keff_history_list`` and ``dominance_ratio()`` are **retired from the
carrier**: :math:`\lambda` is ``sol.outcome.lam`` (and ``outcome.keff``,
which exists only under the k map — under
:data:`~orpheus.numerics.posing.ALPHA_MAP` the same field is
:math:`\alpha` and calling it ``keff`` would be a unit error), the
trajectory is ``outcome.trajectory``, the ratio
``outcome.dominance_ratio()``.  On a ``Solution[SourceOutcome]`` there is
no ``keff`` **attribute at all** — an ``AttributeError``, not a ``None``.

One state, and the flux members READ it
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``outcome.state`` is the returned iterate WHOLE — the one-system coupled
field on a seedless mesh, the two-system one on a carrying (ray-bearing)
mesh — and :attr:`~orpheus.sn.solution.SolutionBase.angular_flux`,
:attr:`~orpheus.sn.solution.SolutionBase.boundary_flux`,
:attr:`~orpheus.sn.solution.SolutionBase.radial_characteristic` and
:attr:`~orpheus.sn.solution.SolutionBase.scalar_flux` are **readers** of
it under their historical names.  Storing the state whole is what lets
three separate guards retire, and each retirement is the same move —
*make the structure say it*:

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - the guard that retired
     - what says it now
   * - a hand-written biconditional asserting the ray member's presence
       matches the Problem's ``R12a`` predicate
     - the state's **ARITY**.  ``radial_characteristic`` is ``None``
       exactly when the state has one system; there is nothing to keep
       in sync
   * - a marginal-axes check keeping a stored ``scalar_flux`` honest
       against the angular one
     - one quantity, one representation.  :math:`\phi = \int\psi\,d\Omega`
       of the state's cell-average moment, cached once per Solution
   * - a space-content check on *each* stored flux field
     - the **state-on-domain law**, once: the state lives on the
       Problem's coupled space, which is the space the recorded question
       is posed on.  A cross-Problem pairing is refused; a same-hub
       cross-KIND ``replace`` is a legal *different solve*, refused by the
       type instead of by a guard

⭐ **A convention divergence closed on the way.**  ``angular_flux`` used
to mean two different things: the eigen and adjoint tails stripped a
multi-moment (LD) closure's :math:`\hat\phi` slopes to the cell average,
while the fixed-source arms kept the whole trailing moment axis.  One
member, two conventions, decided by which entry you called.  Storing the
state whole makes the arm's own convention the *only* convention, and the
cell-average reduction moves to where the scheme lives — the hub's
:meth:`~orpheus.sn.problem.SNProblem.cell_average_moment`, which
the derived :math:`\phi` calls and which is the identity for DD/Step.

.. warning::

   The eigen path's :math:`\phi` is therefore :math:`\int\psi\,d\Omega`
   of the **returned** :math:`\psi` — the one-step-polished,
   section-applied state (:ref:`sn-finalize-one-step`) — and no longer
   the power iteration's own converged scalar.  The two agree to the
   outer residual: ``[M]`` 2026-09-17 against the #448 artefacts as frozen
   BEFORE the step (re-baselined at U2e — the section above),
   ``max rel`` :math:`2.19\times10^{-11}` / :math:`2.94\times10^{-11}` on
   ``cart2d_L0`` / ``slab_vac_L0`` with :math:`k` bit-identical, against
   those gates' :math:`10^{-8}` band (457× and 340× of headroom).  Do
   **not** read that as a floor —
   it is the outer residual, and it falls with the tolerances.

One mint, and the section applied there
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All five public entries and both roles now construct through
:func:`~orpheus.sn.solver._package_solution`.  ⛔ Until step 3 that tail
served three of five while the two fixed-source arms built
``Solution(...)`` inline to keep their DG slope structure — which is
exactly why the convention divergence above was possible, and why no
cross-cutting change at the tail could reach every entry.  Retired with
it: ``_package_adjoint_solution``, ``_exit_balance_defect``,
``_cell_average_angular`` and ``_average_moment_scalar``.

The eigen entries apply the recorded **section** at that mint, so
``gauge.functional(state) == target`` is a LAW of the returned answer
rather than something a consumer must re-check.  It is not a formality:
the driver rescales the SCALAR iterate each outer, and the returned
:math:`\psi` is polished one step past it, so it lands off the section by
the iteration's own residual.  ``[M]`` 2026-09-17, two-region 2-group
slab, ``gauss_legendre(8)``, 8 cells, reading :math:`n(\psi)` **before**
the section: :math:`3.93\times10^{-8}` off at ``solve_sn``'s defaults
(``keff_tol`` 1e-7 / ``flux_tol`` 1e-6 / ``inner_tol`` 1e-8) and
:math:`1.07\times10^{-10}` off at the finalize gates' tolerances
(1e-10 / 1e-9 / 1e-11) — three orders down for three orders of tolerance,
and ``1.0`` to one ulp afterwards.

⚠ **The forward and adjoint eigen sections are DIFFERENT functionals, and
that is now recorded rather than implied.**  The forward records SN's
production rate **including** the :math:`(n,2n)` emission, as a functional
of the STATE; the adjoint records ``KEigenvalue``'s fission-only member.
Both are truthfully described as "unit production rate", and on a
:math:`\Sigma_2`-carrying deck they differ by rel
:math:`1.2\times10^{-1}` — see :ref:`the-outcome-sn-realization` for why
the gauge stores the callable rather than a label.

.. _sn-consuming-the-frame:

Consuming the frame in SN
=========================

Spatial homogenisation and energy condensation are **discrete-frame
projections** — the Petrov-Galerkin coefficient extraction
:math:`G^{-1}M` of a flux- (or spectrum-) weighted frame. All of that
theory — rate preservation, the source-group / sink-sum matrix rules, the
metric-fold-vs-bilinear adjoint argument, fractional-overlap re-binning,
the condense/homogenize asymmetry law, and the verification gates — is
the frame page's headline **Petrov-Galerkin** consumer; see
:ref:`sn-spatial-homogenization` and :ref:`sn-energy-condensation`
(:doc:`/theory/foundations/frame`). This section keeps only the **SN-layer orchestration**:
how the SN :class:`~orpheus.sn.solution.Solution` drives that machinery
from a converged flux.

Homogenisation: the solve → homogenize → re-solve loop
------------------------------------------------------

:meth:`Solution.homogenize <orpheus.sn.solution.Solution.homogenize>`
takes a coarse mesh (:class:`~orpheus.geometry.mesh.Mesh1D` or
:class:`~orpheus.geometry.mesh.Mesh2D`) and returns a
:class:`~orpheus.transport.mesh.material_mesh.MaterialMesh` — the coarse
geometry already carrying one freshly-homogenised effective
:class:`~orpheus.data.macro_xs.mixture.Mixture` per coarse cell. The SN
:class:`~orpheus.sn.solution.Solution` owns the converged flux, so the SN
layer is what builds the flux-weighted **test** basis the frame consumes;
the frame itself, and the rate-preservation theory that *forces* the flux
weighting (rather than a plain volume average), live in
:ref:`sn-spatial-homogenization` (:doc:`/theory/foundations/frame`). The returned
``MaterialMesh`` is re-promoted to a solvable phase space by
:meth:`SNProblem.from_material_mesh
<orpheus.sn.problem.SNProblem.from_material_mesh>`, closing the
**solve → homogenize → re-solve** loop. The return type is
**mesh-coupled** (geometry and materials born together) — the space half
of the condense/homogenize asymmetry law
(:ref:`sn-condense-homogenize-asymmetry`, :doc:`/theory/foundations/frame`).

Condensation: per-material representative spectra
-------------------------------------------------

:meth:`Solution.condense <orpheus.sn.solution.Solution.condense>` is the
SN-layer orchestration of energy condensation. It condenses **each
material with its own representative spectrum** — the flux·volume-weighted
flux over the cells where the material appears:

.. math::
   :label: energy-condensation-representative-spectrum

   \varphi^{(m)}_g \;=\;
   \sum_{i:\,\mathrm{mat}(i)=m} V_i\,\phi_{i,g},

.. (vv-status rationale) Representational identity: the per-material
   representative spectrum used as the condense test weight — the
   flux·volume-weighted flux over the material's cells (mirrors how
   ``homogenize`` derives its flux weight). A definition consumed by
   :meth:`Mixture.condense`; the end-to-end rate preservation it feeds is
   the L1 gate, not a separate claim.
.. vv-status: energy-condensation-representative-spectrum documented

used as the test weight in :meth:`Mixture.condense
<orpheus.data.macro_xs.mixture.Mixture.condense>` — the data-layer
collapse verb, whose spectrum-weighted-collapse theory is
:ref:`sn-energy-condensation` (:doc:`/theory/foundations/frame`) — mirroring how
:meth:`Solution.homogenize` derives its flux weight from the same solved
flux. The result is a **portable** ``dict[int, Mixture]`` keyed by
material id — few-group cross sections carrying the coarse ``eg``, not
bound to any mesh (the **mesh-decoupled** half of the asymmetry law,
:ref:`sn-condense-homogenize-asymmetry`, :doc:`/theory/foundations/frame`). A material with
no flux in a fine group contributes zero weight there; the condense
frame's Moore–Penrose Gram handles any empty coarse group.
