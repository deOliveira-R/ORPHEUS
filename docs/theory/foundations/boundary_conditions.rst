.. _theory-boundary-conditions:

============================================
Boundary Conditions — Trace-Law Architecture
============================================

.. contents:: Contents
   :local:
   :depth: 3


Key Facts
=========

**Read this before touching anything in** :mod:`orpheus.geometry.boundary`
**or any** ``BoundaryRealizer``.

- A boundary condition is a **method-agnostic affine map** on the
  transport equation's boundary trace:
  :math:`\gamma_- \psi = R\,G\,\gamma_+ \psi + q`, where
  :math:`\gamma_\pm` are the inflow / outflow trace operators,
  :math:`G` is a geometric map (permutation, pushforward, angular
  average), :math:`R` is a response amplitude (:term:`albedo`, sub-Markov
  kernel), and :math:`q` is an optional prescribed inflow source.
  See :eq:`affine-bc-form`.
- The architecture has **three concrete layers**, connected by the
  kind-keyed law registry (#290 P7b dissolved the Wave-5 realizer
  registry — realizers are owned by their method-meshes):

  +-------+-----------------------+-----------------------------------------------------+
  | Layer | What                  | Where                                               |
  +=======+=======================+=====================================================+
  | 1     | Trace structure       | :mod:`orpheus.numerics.spaces.angular_trace_space`  |
  |       | (Γ\_-, Γ\_+ + mask)   | (all Mesh1D coord systems + 2-D Cartesian;          |
  |       |                       | 2-D cylindrical Mesh2D deferred)                    |
  +-------+-----------------------+-----------------------------------------------------+
  | 2     | Boundary law          | :mod:`orpheus.geometry.boundary` (ABC +             |
  |       | (method-agnostic)     | 7 concrete laws, kind-keyed law registry)           |
  +-------+-----------------------+-----------------------------------------------------+
  | 3     | Method realizer       | per-method packages (SN + diffusion,                |
  |       | (per-method strategy) | #290 P3), each owned by its method-mesh             |
  |       |                       | via the ``TransportMethod`` hook (P7b)              |
  +-------+-----------------------+-----------------------------------------------------+

- Rank-N boundary conditions (Marshak, partial-current mixes) are
  expressed via a **descriptor-tree algebra** on the unrealised laws
  themselves. The :class:`~orpheus.geometry.boundary.BoundaryTraceLaw`
  algebra dunders (``+``, ``-``, ``*``, ``/``, unary ``-``) return
  :class:`~orpheus.geometry.boundary.LawSum` /
  :class:`~orpheus.geometry.boundary.LawScaled` nodes — a closed
  algebra over ``BoundaryTraceLaw | LawSum | LawScaled``. The tree is
  a **pure descriptor** with no ``apply`` method; the
  :func:`~orpheus.geometry.boundary.realize_recursively` type
  transformer (method-blind, in ``geometry/boundary/`` since #290
  P7b; the leaf realizer is a required argument) walks it once per
  face and emits an operator tree of
  :class:`~orpheus.numerics.operator.OperatorSum` /
  :class:`~orpheus.numerics.operator.ScaledOperator` composers around
  realised 1-arg leaves. See :ref:`bc-trace-law-descriptor-model` and
  :ref:`bc-rank-n-algebra`. There is no dedicated
  ``MixedBoundaryOperator`` class (retired Wave 11); there is also no
  ``apply`` method on the raw law (retired Issue #186, B3 + β2).
- TWO functional realizers exist:
  :class:`~orpheus.sn.boundary.realizer.SNBoundaryRealizer` and
  :class:`~orpheus.diffusion.boundary_realizer.DiffusionBoundaryRealizer`
  (#290 P3 — every diffusion law collapses to the albedo-family
  scalar :math:`\mathcal{A}` in :math:`J^- = \mathcal{A} J^+`). Each
  is OWNED by its method-mesh: ``realize_boundary_law`` — the
  per-method arm of the
  :class:`~orpheus.transport.method.TransportMethod` Protocol —
  instantiates it directly, and the shared
  :func:`~orpheus.transport.method.resolve_boundary_conditions` body
  drives the per-face resolution for every method. The Wave-5
  ``BoundaryRealizerRegistry`` + the MoC/MC/CP
  ``NotImplementedError`` stub realizers were **dissolved at #290
  P7b**: no consumer ever resolved a realizer by name (you hold the
  method-mesh → you have its realizer), and the string indirection
  carried registration-timing hazards for zero payoff. A future
  MoC / MC / CP modernization mints its method-mesh + realizer pair
  directly — no central registration step.
- The :attr:`creates_sweep_cycle` ``ClassVar`` flag on each
  :class:`~orpheus.geometry.boundary.BoundaryTraceLaw` subclass
  signals to the SN :term:`sweep` planner (§15A.2) which boundary types
  introduce cycles in the directed cell-visit graph. ``True`` on
  :class:`~orpheus.geometry.boundary.ReflectiveBoundary` and
  :class:`~orpheus.geometry.boundary.PeriodicBoundary`; ``False``
  on all other laws.
- The eight typed errors :class:`~orpheus.geometry.boundary.IncomingOutgoingTraceClassificationError`
  through :class:`~orpheus.geometry.boundary.BoundarySourceNotOnIncomingTraceError`
  (ERR-040..ERR-047 in the V&V error catalog at
  ``.claude/skills/vv-principles/error_catalog.md``) replace the
  pre-refactor generic :class:`ValueError` raises; every one is
  pinned by a ``@pytest.mark.catches("ERR-NNN")`` decorator on the
  test that fires it.
- **Vacuum semantic correction (§16A.5).** ``VacuumInflow`` realises
  to :class:`~orpheus.numerics.operator.IncomingOrdinateMaskTensor`,
  which zeroes **only the inflow ordinates** and preserves the
  outflow trace. The pre-refactor ``VacuumBoundaryOperator.apply``
  used ``np.zeros_like(psi_out)`` and zeroed everything; the §16A.10
  inflow-only mask is the trace-correct representation.
  **Post Issue #186 (2026-05-11)** the realizer-routed
  inflow-only mask is the **sole** path to vacuum action — the
  zeros-all body has been deleted along with every other
  :meth:`apply` method on every concrete BC. The §16A.5
  two-paths-divergence is therefore eliminated by design (no
  second path remains). The realizer-routed mask is uniform across
  every supported mesh (1-D Cartesian / spherical / cylindrical +
  2-D Cartesian) since Issue #188 lifted the curvilinear deferral on
  the boundary trace space (then named ``InflowTraceSpace``; since
  #205 / #201 the unified
  :class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace`).
- **The realized boundary law is a first-class sibling operator**
  :math:`B` **in the SN algebra (Wave O steps O.4a.2 + O.4b, Issue
  #208).** For **every** SN geometry (1-D slab / sphere / cylinder and
  2-D Cartesian), the realized per-face law is assembled into the
  whole-trace
  :class:`~orpheus.sn.operators.boundary.SNBoundaryOperator` — the
  :math:`A_{ss}` boundary block of the canonical loss operator
  :math:`(L_{\rm full} + C - S - F - B)`. The reflection
  :math:`\psi.\text{inflow} = B\,\psi.\text{outflow}` is **no longer
  re-applied inside the streaming sweep** for any geometry (O.4a.2
  made the 1-D sweep bare; O.4b made the 2-D wavefront sweep + matvec
  bare); it is delivered as the off-diagonal :math:`-B` source term and
  the outer Krylov / SI loop drives the boundary consistency residual
  to zero. The full block-matrix derivation and design rationale live
  at :ref:`bc-extraction`.

.. admonition:: V&V status

   This page is **L4-informational** with respect to correctness.
   The architecture documented here is structural — it makes the
   code understandable and composable but does not by itself verify
   any equation. The verification load is carried by:

   - L0 foundation tests on individual primitives
     (:mod:`tests.numerics`,
     :mod:`tests.geometry.test_boundary_trace_law`,
     :mod:`tests.geometry.test_bc_errors`).
   - L1 realiser-output snapshot tests (the Wave-6 harness at
     :mod:`tests.geometry.test_bc_equivalence_snapshot`, with the
     legacy halves dropped post Issue #186 / C-B3.7). The
     surviving ``test_realizer_*`` halves pin the realised-operator
     output against committed ``.npz`` snapshots at ``nulp ≤ 4``
     for non-vacuum BCs (intentional semantic capture of the
     §16A.5 inflow-only mask for vacuum).
   - L1 descriptor-tree algebra tests
     (:mod:`tests.geometry.test_law_composition`) pinning the
     :class:`LawSum` / :class:`LawScaled` closed-algebra contract
     (foundation + L1 coverage).
   - L1 universal-invariant tests
     (:mod:`tests.geometry.test_bc_universal_invariants`) that fire
     ERR-043 / ERR-044 / ERR-046 under fault-injection.

   No equation on this page makes a claim that requires a closed-form
   or MMS reference; all equations are **definitional** or
   **structural-architecture** statements drawn from Grand Report v3
   §16, §16A, §16A.10.


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
   PrescribedInflow) are the tested forms. Not a solver claim.
.. vv-status: affine-bc-form documented

where:

* :math:`\gamma_\pm` are the **trace operators** that restrict the
  angular flux :math:`\psi(\mathbf{r}, \Omega)` from the volumetric
  function space to the inflow / outflow boundary trace spaces
  :math:`\Gamma_\pm` (see :ref:`bc-trace-structure` below for the
  formal definition).
* :math:`G : \Gamma_+ \to \Gamma_+` is the **geometric map** — a
  measure-preserving permutation, pushforward, spatial wrap-around,
  or hemispheric cosine-weighted average. It carries pure geometry
  (it changes nothing about the physical interaction at the
  boundary; it just relabels the angular fluxes that meet there).
* :math:`R : \Gamma_+ \to \Gamma_-` is the **response kernel** — a
  scalar amplitude in :math:`[0, 1]` for the standard sub-Markov BCs
  (albedo, white, partial-current) or a full angular kernel in
  general weak-form BCs (deferred; see the
  :class:`~orpheus.geometry.boundary.BoundaryError` catalog and
  Issue #175 close-out follow-ups).
* :math:`q \in \Gamma_-` is the **prescribed inflow source** — a
  vector-valued quantity on :math:`\Gamma_-` only. The empty case
  :math:`q \equiv 0` is the homogeneous BC; the inhomogeneous case
  :math:`q \neq 0` is the rank-0 affine BC
  :class:`~orpheus.geometry.boundary.PrescribedInflow`.

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
   form an explicit name; ORPHEUS does because two distinct rank-0
   cases (:class:`~orpheus.geometry.boundary.VacuumInflow` with
   :math:`R = G = q = 0` and
   :class:`~orpheus.geometry.boundary.PrescribedInflow` with
   :math:`R = G = 0` but :math:`q \neq 0`) need a single uniform
   contract.
3. **The three operators are first-class.** The
   :class:`~orpheus.geometry.boundary.BoundaryTraceLaw` ABC exposes
   :attr:`~orpheus.geometry.boundary.BoundaryTraceLaw.geometry_map`,
   :attr:`~orpheus.geometry.boundary.BoundaryTraceLaw.response_kernel`,
   and :attr:`~orpheus.geometry.boundary.BoundaryTraceLaw.source`
   as Python properties on every concrete subclass. The properties
   default to ``None / None / NoSource()``; concrete laws override
   when applicable. The split lets cross-method realizers introspect
   the law's geometric and response components separately — the SN
   realizer dispatches on the law's class today, but a future
   weak-form realizer might dispatch on the geometry / response /
   source independently.

.. note::

   **SN apply matvec honours the affine BC contract (Issue #168
   Phase C, 2026-05-12).** The then-production per-geometry matvecs
   (``transport_operator_matvec_spherical`` and ``_cylindrical`` —
   since deleted in the typed-field campaign (#197), their successor
   ``_transport_operator_matvec_unified`` in turn retired at the walk
   unification (#280 campaigns)) were rewritten as one sweep
   iteration semantically: the BC trace law is applied **at least
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
   (Phase D Carlson context + Phase C trace law).


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
  is constrained by the law.
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

The signed-projection table is what the
:class:`~orpheus.sn.boundary.realizer.SNBoundaryRealizer` reads to build
the sparse vacuum-mask operator
(:class:`~orpheus.numerics.operator.IncomingOrdinateMaskTensor`) that
zeros precisely the inflow ordinates at one face.


.. _bc-realizer-layer:

Layer 3 — the method realizer
-----------------------------

A single :class:`~orpheus.geometry.boundary.BoundaryTraceLaw`
describes the physics at the boundary but is **not** by itself
ready for consumption by a transport sweep. The conversion from
method-agnostic law to method-specific
:class:`~orpheus.numerics.operator.LinearOperator` is the job of a
:class:`~orpheus.geometry.boundary.BoundaryRealizer`. Each transport
method that has adopted the unified BC architecture (SN; diffusion
since #290 P3) ships one realizer class, owned by its method-mesh:
the mesh's ``realize_boundary_law`` arm — the per-method hook of the
:class:`~orpheus.transport.method.TransportMethod` Protocol (#290
P7b) — instantiates it directly.

The realizer's :meth:`realize` method takes the law plus a
**method space** — a method-specific container holding the
:term:`quadrature`, mesh, trace masks, and any other discretization
metadata the realizer needs — and returns a 1-arg
:class:`~orpheus.numerics.operator.LinearOperator` whose
:meth:`apply` carries the method-specific realization of the
affine BC :eq:`affine-bc-form`.

Why this third layer? Because the same affine law is realized by
**structurally different** linear operators in each transport
method:

* SN realizes vacuum as a sparse per-ordinate
  :class:`~orpheus.numerics.operator.IncomingOrdinateMaskTensor` on
  the inflow ordinates of the affected face.
* MoC realizes vacuum by zeroing the entering boundary fluxes of
  every track that intersects the face.
* MC realizes vacuum by killing particle histories at the face.
* CP realizes vacuum as zero rows in the boundary-to-region
  coupling matrix.

Splitting the realizer out of the law makes each piece independently
testable and gives every method a single bolt-in point for its BC
treatment — see :ref:`bc-cross-method-stubs` for how a future method
adopts the architecture (and for the Wave-5 stub scaffolding that was
dissolved at #290 P7b).


.. _bc-extraction:

Boundary-condition extraction — :math:`B` as a sibling operator
===============================================================

The composite metric-correct G-adjoint that closes this extraction
narrative is documented at :ref:`g-adjoint` in
:doc:`/theory/foundations/operator_adjoint`.

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
prediction (:ref:`tensor-network-decomposition`): "Wave O typing must accept
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
  :math:`A_{ss}` (:math:`R` the realized per-face law,
  :ref:`bc-law-layer`) sits **strictly upper** on the same ordering —
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
       realized per-face law :math:`\psi.\text{inflow} =
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

The honest :math:`L+C-S-B` driver via variadic couplings
--------------------------------------------------------

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
(:class:`~orpheus.sn.operators.streaming.StreamingCollisionOperator`, ``.solve`` = the WDD
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
**same** predicate the representation dispatch
(:func:`~orpheus.sn.loss_representation.default_for`, via each
representation's ``supports``) reads to select the 1-D scan
(:class:`~orpheus.sn.loss_representation.CumprodScan`) vs the 2-D
wavefront body, and the **same**
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
       (``tests/sn/sweep/curvilinear/test_streaming_equilibrium_curvilinear.py``)
     - Analytical infinite-medium balance
       :math:`\phi = q/\Sigma_t` (closed-form)
     - ``source_iteration`` AND ``krylov``
   * - Reflective :math:`k_\infty` homogeneous
       (``tests/sn/verification/analytical/test_kinf_homogeneous.py``)
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

Operator-output role typing — :math:`A\psi` is a source/sink
------------------------------------------------------------

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
  V_{\rm trace}` (:ref:`g-adjoint`). This is exactly why
  :math:`B` stays trace-typed as a separate gain
  (:ref:`bc-extraction-variadic-driver`): a bulk-folded :math:`B` would
  erase the trace metric the adjoint needs.
* **Gate-1.3** (the O.2 adjoint verification gate) **landed in R5** —
  the dense-probe oracle + the L11 wrong-metric control
  (:ref:`g-adjoint`).

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
:math:`\Omega \cdot \hat n`:

.. math::
   :label: trace-sign-predicate

   \Gamma_- \;=\; \{(\mathbf{r}, \Omega) \in \partial\Omega \times S^d
                  : \Omega \cdot \hat n(\mathbf{r}) < 0\},
   \qquad
   \Gamma_+ \;=\; \{(\mathbf{r}, \Omega) \in \partial\Omega \times S^d
                  : \Omega \cdot \hat n(\mathbf{r}) > 0\}.

.. (vv-status rationale) Notation definition: the continuous inflow / outflow
   trace half-spaces Γ_± by the sign of Ω·n. Its discrete realisation
   :eq:`inflow-mask-discrete` is the tested form
   (``tests/numerics/test_angular_trace_space.py`` selector gates). A
   definitional predicate, not a solver claim.
.. vv-status: trace-sign-predicate documented

Points with :math:`\Omega \cdot \hat n = 0` are **tangential** —
they belong to neither half. For axis-aligned ordinates on
axis-aligned faces these arise exactly (no round-off) for face
normals perpendicular to the ordinate's direction cosine; for
general curvilinear faces or generic ordinates they arise only
at a measure-zero subset that the discrete representation
identifies via a small tolerance (``_TANGENTIAL_EPS = 1e-12``).

In the discrete setting, the spatial boundary is a union of finite
faces :math:`\{f_1, \ldots, f_F\}` and the angular variable is a
finite ordinate set :math:`\{\Omega_n : n = 1, \ldots, N\}`. The
sign predicate :eq:`trace-sign-predicate` then collapses to a
**per-face boolean mask** of shape :math:`(F, N)`:

.. math::
   :label: inflow-mask-discrete

   \mathrm{inflow\_mask}[f, n]
   \;=\; \bigl(\Omega_n \cdot \hat n_f < -\epsilon\bigr),
   \qquad
   \mathrm{outflow\_mask}[f, n]
   \;=\; \bigl(\Omega_n \cdot \hat n_f > +\epsilon\bigr).

This mask is the discrete realization of :math:`\Gamma_\pm`. It is
the load-bearing primitive that downstream consumers need:

* The SN realizer's vacuum branch reads
  ``inflow_mask[f]`` for the specific face :math:`f` and converts
  it to an integer array of ordinate indices via
  :meth:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.inflow_indices_for_face`.
  Those indices are the constructor argument to
  :class:`~orpheus.numerics.operator.IncomingOrdinateMaskTensor`.
* The universal invariant
  :meth:`~orpheus.geometry.boundary.BoundaryTraceLaw.assert_source_lives_on_incoming_trace`
  uses the inflow mask to validate that a
  :class:`~orpheus.geometry.boundary.InflowSourceSpec` has no nonzero
  entries on outflow ordinates (per
  :class:`~orpheus.geometry.boundary.BoundarySourceNotOnIncomingTraceError`,
  ERR-047).
* The SN curvilinear sweep (1-D spherical / cylindrical) consumes
  the same realizer-routed mask as the Cartesian path — Issue #188
  (C188.1+C188.2 in :mod:`orpheus.numerics.spaces.angular_trace_space`, C188.3 in
  :mod:`orpheus.sn.mesh.augmented_mesh`) lifted the curvilinear deferral and
  Issue #176 then dropped the legacy 2-arg shim that existed only
  to bridge that deferral.

The class :class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace`
carries the per-face :math:`\Omega\cdot\hat n` masks as
``Optional[np.ndarray]`` fields excluded from
:meth:`__eq__` and :meth:`__hash__` — preserving the
:class:`~orpheus.numerics.space.FunctionSpace` identity convention
``(name, shape)``. Construction goes through the classmethod factory
:meth:`AngularTraceSpace.from_quadrature_and_layout
<orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.from_quadrature_and_layout>`;
the bare dataclass constructor is reserved for trace spaces whose mask
will be populated later (or never).

**Geometry-blind, layout-driven (post Issue #225 / C5.3).** The factory
is **geometry-blind**: it takes only the angular quadrature and a
:class:`~orpheus.numerics.face_layout.FaceLayout` (canonically
:attr:`SNMesh.boundary_face_layout <orpheus.sn.mesh.augmented_mesh.SNMesh>`), and
reads every datum from those two — the layout's ``"{axis}{min|max}"``
face names imply axis-aligned outward normals, so the
:math:`\Omega\cdot\hat n` row for an axis-:math:`a` face is
:math:`\pm\mu_a` (sign from ``min`` / ``max``). It works for every
constructible :class:`~orpheus.geometry.mesh.Mesh1D` coord system
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
such a :class:`~orpheus.geometry.mesh.Mesh2D` cannot become an
:class:`SNMesh` (no 2-D cylindrical SN sweep exists); the refusal lives
at the :class:`SNMesh` construction surface, not in the trace factory.


.. _bc-law-layer:

Boundary law (``BoundaryTraceLaw`` ABC + concretes)
===================================================

The base class :class:`~orpheus.geometry.boundary.BoundaryTraceLaw`
is an ``abc.ABC`` that combines
:class:`~orpheus.numerics.operator.LinearOperator` (for
operator-algebra dunders like ``+``, ``*``, ``@``) and
:class:`~orpheus.numerics.registry.RegistryMixin` (so each concrete
subclass self-registers under its ``key=`` class-creation kwarg).
The ABC ships:

1. Three first-class properties carrying the
   :eq:`affine-bc-form` operators: ``geometry_map``,
   ``response_kernel``, ``source`` (defaulting to ``None``, ``None``,
   :class:`~orpheus.geometry.boundary.NoSource`).
2. The :attr:`~orpheus.geometry.boundary.BoundaryTraceLaw.creates_sweep_cycle`
   ``ClassVar`` boolean signal — see :ref:`bc-sweep-cycle`.
3. Five **universal** ``assert_*`` invariants (no-op defaults;
   concrete laws override the relevant ones) and four **specific**
   invariants on the BCs that need them — see
   :ref:`bc-universal-invariants`.
4. A :meth:`realize` hook that raises with guidance (route through a
   method realizer) — see :ref:`bc-realizer-layer-detail`.
5. **No ``apply`` method at all** (Issue #186 / B3 + β2,
   2026-05-11). The descriptor model that survived the C176.3
   Option A interim was retired in favour of a pure-descriptor
   contract: :class:`BoundaryTraceLaw` is no longer a
   :class:`~orpheus.numerics.operator.LinearOperator` subclass,
   no concrete law carries ``apply`` / ``apply_transpose``
   methods, and none reports the operator predicates
   ``is_invertible`` / ``is_adjointable``.
   The §16A.3 three-layer split (descriptor / realizer / operator)
   is now enforced by the **type system**, not by convention —
   ``law.apply(psi)`` on a raw law is an ``AttributeError`` at
   runtime and a static-type error at the linter level. The
   :class:`~orpheus.sn.boundary.realizer.SNBoundaryRealizer` is the
   sole bridge from descriptor to callable; see
   :ref:`bc-trace-law-descriptor-model` for the design rationale
   and the predecessor approaches that were tried and rejected.

The six concrete laws ship under :mod:`orpheus.geometry.boundary`,
one per submodule. The Grand Report v3 vocabulary is used verbatim
for the class names; pre-refactor names are retained as deprecated
aliases in the package ``__init__`` (see :ref:`bc-naming-audit`).

.. list-table:: Concrete ``BoundaryTraceLaw`` subclasses
   :header-rows: 1
   :widths: 18 16 13 13 13 27

   * - Class
     - Registry key
     - :math:`G_\alpha`
     - :math:`R_\alpha`
     - :math:`q`
     - Sweep-cycle flag
   * - :class:`~orpheus.geometry.boundary.VacuumInflow`
     - ``"vacuum"``
     - rank-0 (none)
     - 0
     - 0
     - ``False``
   * - :class:`~orpheus.geometry.boundary.ReflectiveBoundary`
     - ``"reflective"``
     - axis-reflection permutation
     - albedo
     - 0
     - **``True``** (signals §15A.2)
   * - :class:`~orpheus.geometry.boundary.WhiteBoundary`
     - ``"white"``
     - cosine-weighted hemispheric average (Lambertian)
     - albedo
     - 0
     - ``False``
   * - :class:`~orpheus.geometry.boundary.PeriodicBoundary`
     - ``"periodic"``
     - spatial pushforward (caller-supplied)
     - 1
     - 0
     - **``True``** (signals §15A.2)
   * - :class:`~orpheus.geometry.boundary.AlbedoBoundary`
     - ``"albedo"``
     - identity in angle
     - albedo
     - 0
     - ``False``
   * - :class:`~orpheus.geometry.boundary.PrescribedInflow`
     - ``"prescribed_inflow"``
     - 0
     - 0
     - :math:`q \in \Gamma_-`
     - ``False``

Per-law rank census (each kernel read off the :math:`G_\alpha`
column above): vacuum and prescribed inflow carry no response
kernel (rank-0). White is **rank-1 in angle** — one isotropic
re-entry mode, fed by the cosine-weighted (Lambertian) outflow
average, written to every inflow ordinate. Reflective and albedo
are **angular permutations** (axis-reflection and identity
respectively, scaled by the albedo factor) — rank :math:`N/2` per
slab face: structured, trace-sized, NOT rank-1. Periodic is a
**spatial pushforward** pairing opposite faces. Marshak /
partial-current boundaries are rank-N via the
**descriptor-tree algebra** on the unrealised laws (:class:`LawSum`
/ :class:`LawScaled` over :class:`BoundaryTraceLaw` leaves) —
realised once per face by
:func:`~orpheus.geometry.boundary.realize_recursively`; see
:ref:`bc-rank-n-algebra` below. Every one of these kernels is
**trace-sized** — tiny against the bulk dimension — which is
exactly the low-rank/structured shape a Woodbury boundary closure
exploits (:ref:`the scoped SMW statement <smw-low-rank-exception>`,
Issue #300).


.. _bc-naming-audit:

Naming audit: pre-refactor vs Grand Report v3 vocabulary
--------------------------------------------------------

Wave 7 of the refactor renamed every concrete BC to match the
Grand Report v3 vocabulary verbatim. During the deprecation window
the pre-refactor names were re-exported as deprecated aliases from
:mod:`orpheus.geometry.boundary.__init__` so existing import sites
kept working unchanged. Those aliases were **retired in Wave O step
O.4a.1** once every code and test consumer had migrated; the
canonical names in the middle column are now the sole importable
symbols. The table is retained as the historical naming index for
readers tracing pre-Wave-O commits:

.. list-table:: Wave 7 BC renames (pre-refactor name → canonical name)
   :header-rows: 1
   :widths: 35 35 30

   * - Pre-refactor name (retired Wave O O.4a.1)
     - Canonical name
     - Why renamed
   * - ``VacuumBoundaryOperator``
     - :class:`~orpheus.geometry.boundary.VacuumInflow`
     - Emphasizes "inflow set to zero", not "operator that vacuums";
       distinguishes from the rank-N case
       :class:`PrescribedInflow` which also writes only the inflow
       trace.
   * - ``SpecularBoundaryOperator``
     - :class:`~orpheus.geometry.boundary.ReflectiveBoundary`
     - "Specular" is one specific axis-aligned reflection;
       "Reflective" is the family name that the Grand Report uses.
       A future ``SymmetryBoundary`` (deferred) will share the
       reflective-family base but apply on non-physical octant
       boundaries.
   * - ``WhiteBoundaryOperator``
     - :class:`~orpheus.geometry.boundary.WhiteBoundary`
     - Drops the redundant "Operator" suffix that pre-dated the
       law / realizer split. The law is no longer "the operator";
       it's the abstract physical statement that gets *realized*
       to an operator.
   * - ``PeriodicBoundaryOperator``
     - :class:`~orpheus.geometry.boundary.PeriodicBoundary`
     - Same rationale: pre-refactor "Operator" suffix is
       structurally misleading.
   * - ``AlbedoBoundaryOperator``
     - :class:`~orpheus.geometry.boundary.AlbedoBoundary`
     - Same rationale.
   * - ``MixedBoundaryOperator``
     - **retired Wave 11** — see :ref:`bc-rank-n-algebra`.
     - Replaced by the descriptor-tree algebra
       (:class:`LawSum` / :class:`LawScaled` over
       :class:`BoundaryTraceLaw` leaves; Issue #186 / B3 + β2);
       the dedicated class added no value over the inherited
       algebra dunders.


.. _bc-sweep-cycle:

The ``creates_sweep_cycle`` signal
----------------------------------

The SN sweep is a topological sort of the directed cell-visit graph
where edges are oriented by :math:`\mathrm{sign}(\Omega_n \cdot
\hat n_f)`. For most BCs the boundary is the *root* of the sort —
inflow values come from the BC, get propagated through the cells,
and exit as outflow values that the BC consumes but doesn't feed
back. For two BC families this is no longer true:

* **Reflective.** The outflow flux at a face is mapped to an inflow
  flux at the **same** face (under the reflection permutation). The
  cell visits that produce the outflow are predecessors of the cell
  visits that consume the reflected inflow — but those are *the
  same cells*. The sweep DAG acquires a cycle and a self-consistent
  fixed-point iteration must converge the reflective inflow against
  the outflow.
* **Periodic.** Same situation, but the cycle spans two different
  faces (outflow at face A maps to inflow at face B, and vice
  versa).

The :attr:`~orpheus.geometry.boundary.BoundaryTraceLaw.creates_sweep_cycle`
``ClassVar`` is the structural signal that lets the §15A.2
sweep-cycle detector identify these BCs without inspecting the
realization. It is ``True`` on :class:`ReflectiveBoundary` and
:class:`PeriodicBoundary`, ``False`` on every other law.

Vacuum, white, albedo, and prescribed-inflow are all cycle-free:
they consume :math:`\gamma_+ \psi` and produce a function of it
that only depends on :math:`\gamma_+` (or, for prescribed inflow,
nothing at all). No fixed point is needed.


.. _bc-realizer-layer-detail:

Realizer (``BoundaryRealizer`` Protocol + ``SNBoundaryRealizer``)
=================================================================

The :class:`~orpheus.geometry.boundary.BoundaryRealizer` Protocol is
``@runtime_checkable`` and lives at
:mod:`orpheus.geometry.boundary._realizer`. Its contract is one
attribute and one method:

.. code-block:: python

   @runtime_checkable
   class BoundaryRealizer(Protocol):
       method_name: str
       def realize(
           self,
           law: BoundaryTraceLaw,
           method_space: Any,
       ) -> LinearOperator: ...

The Protocol intentionally does *not* prescribe how
:meth:`realize` dispatches over law types — different methods will
have different optimal dispatch strategies. The SN realizer uses
``isinstance`` because the law set is small and stable; a future
realizer that needs runtime extension might use the
:class:`~orpheus.numerics.registry.RegistryMixin` machinery instead.

The Wave 5 SN dispatch table is the documented standard — the §15.2
:math:`G_\alpha` geometric primitives, one per boundary law:

.. _bc-tensor-primitives:

.. list-table:: SN realization map (law → Wave-0 / Wave-1 primitive)
   :header-rows: 1
   :widths: 24 38 38

   * - Law
     - Realized representation (α = 1 fast path)
     - Realized representation (α ∉ {0, 1})
   * - :class:`VacuumInflow`
     - :class:`~orpheus.numerics.operator.IncomingOrdinateMaskTensor`
       with the per-face ``inflow_indices`` from the method space's
       :class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace`
       (selected by
       :meth:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.inflow_indices_for_face`).
     - n/a (vacuum has no α parameter)
   * - :class:`ReflectiveBoundary(axis, α)`
     - bare
       :class:`~orpheus.numerics.operator.PermutationOperator`
       with ``perm = quadrature.reflection_index(axis)``
     - ``ScaledOperator(α, PermutationOperator(...))``
   * - :class:`WhiteBoundary(axis, outward_sign, α)`
     - bare
       :meth:`~orpheus.sn.boundary.angular.AngularAverageOperator.from_quadrature`
       (``from_quadrature(quad, axis, outward_sign)``)
     - ``ScaledOperator(α, AngularAverageOperator(...))``
   * - :class:`AlbedoBoundary(α)` with α=0
     - :class:`~orpheus.numerics.operator.ZeroOperator`
     -
   * - :class:`AlbedoBoundary(α)` with α=1
     - :class:`~orpheus.numerics.operator.IdentityOperator`
     -
   * - :class:`AlbedoBoundary(α)` with α ∉ {0, 1}
     -
     - ``ScaledOperator(α, IdentityOperator())``
   * - :class:`PeriodicBoundary`
     - :class:`~orpheus.numerics.operator.PeriodicWrapOperator`
       (today an angular identity; spatial-pushforward extension
       pending — see "BC: PeriodicWrapOperator spatial-pushforward
       implementation" follow-up, ``module:sn``)
     - n/a (periodic has no α parameter)
   * - :class:`PrescribedInflow(source)`
     - :class:`~orpheus.sn.boundary.angular.IncomingSourceOperator`
       (``IncomingSourceOperator(source)``)
       — :meth:`apply` ignores the outgoing flux and returns
       ``source.evaluate(probe_inflow_trace)``
     - n/a

The α = 1.0 fast paths return the **bare** primitive (no
``ScaledOperator`` wrap). This is load-bearing for bit-identity:
without it, the "perfect reflection" case
:class:`~orpheus.geometry.boundary.ReflectiveBoundary` (pre-refactor
``SpecularBoundaryOperator(axis="x", albedo=1.0)``) would shift by
one ULP under the realizer relative to its pre-refactor
``np.take(psi_out, reflection_index, axis=0)`` body — see the
Wave 6 snapshot harness for the bit-equivalence pin.

The :class:`~orpheus.sn.mesh.method_space.SNMethodSpace` dataclass is the
realizer's second argument. It carries:

* :attr:`~orpheus.sn.mesh.method_space.SNMethodSpace.quadrature` — the
  angular quadrature (mandatory).
* :attr:`~orpheus.sn.mesh.method_space.SNMethodSpace.face` — the
  face-name label (``"xmin"`` … ``"zmax"``) so the vacuum branch can
  look up the right inflow indices. (The pre-C4 ``"left"`` / ``"right"``
  spellings were aliases of ``"xmin"`` / ``"xmax"``; since the C4
  face-name carve — :ref:`bc-face-name-carve` — every face is keyed by
  its canonical ``"{axis}{min|max}"`` name.)
* :attr:`~orpheus.sn.mesh.method_space.SNMethodSpace.inflow_indices` —
  the per-face inflow ordinate indices for the vacuum branch
  (derived from the held trace at :meth:`for_face` time).
* :attr:`~orpheus.sn.mesh.method_space.SNMethodSpace.mesh`,
  :attr:`~orpheus.sn.mesh.method_space.SNMethodSpace.trace` — the
  (optional) spatial mesh and the single unified
  :class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace` for any
  realizer branch that needs more than the per-face index list. Since
  the #205 / #201 unification this is **one** ``trace`` attribute, not a
  separate ``inflow_trace`` / ``outflow_trace`` pair — inflow and
  outflow are selectors over the same trace space (see the
  one-space-two-selectors note above). The ``mesh`` slot is optional
  metadata (C5.3, #225): nothing in the realizer chain reads it (the
  inflow indices come from the trace), and an axis-native ``SNMesh``
  with no legacy mesh adapter passes ``None``.

The :meth:`SNMethodSpace.for_face` factory is the standard
construction site inside ``SNMesh.realize_boundary_law`` (per face,
driven by the shared
:func:`~orpheus.transport.method.resolve_boundary_conditions` body,
#290 P7b); the :meth:`SNMethodSpace.minimal` factory returns a
quadrature-only method space for unit tests that don't need mesh +
face metadata.


.. _bc-dual-registry:

The law registry (and the realizer registry that was dissolved)
===============================================================

ONE registry connects the tag layer to the law layer:

**Law registry** — keyed by ``BC.kind`` string
(``"vacuum"``, ``"reflective"``, ``"white"``, ``"periodic"``,
``"albedo"``, ``"prescribed_inflow"``, ``"zero_flux"``). The
registry IS :attr:`BoundaryTraceLaw.registry` (a class-level dict
maintained by :class:`~orpheus.numerics.registry.RegistryMixin`).
Concrete laws self-register at module import time via the
``key=`` class-creation kwarg:

.. code-block:: python

   class VacuumInflow(BoundaryTraceLaw, key="vacuum"):
       ...

Lookup is :meth:`BoundaryTraceLaw.create("vacuum")` or direct
dictionary access ``BoundaryTraceLaw.registry["vacuum"]``. Each
method-mesh additionally carries its own **admission table**
(``BOUNDARY_OPERATOR_REGISTRY: dict[str, type[BoundaryTraceLaw]]``
— the subset of laws its realizer can honestly realize, e.g.
``zero_flux`` is diffusion-only), and the shared
:func:`~orpheus.transport.method.resolve_boundary_conditions` body
uses THAT table to recover the law class from a mesh-declared
:class:`~orpheus.geometry.mesh.BC` — an unsupported tag refuses at
phase-space construction with the method's supported list.

**The realizer registry was dissolved at #290 P7b.** The Grand
Report §16A.11 design (lines 3252–3257) paired the law registry
with a second, method-name-keyed ``BoundaryRealizerRegistry``
(``"SN"``, ``"MoC"``, …; realizers self-registered via a decorator
at import time). It shipped in Wave 5 and was retired when the
second functional realizer (diffusion, #290 P3) made the real
consumption pattern visible: **no consumer ever resolved a realizer
by method name**. Production holds a method-mesh, and the mesh's
``realize_boundary_law`` arm — the per-method hook of the
:class:`~orpheus.transport.method.TransportMethod` Protocol —
instantiates its own realizer directly; the rank-N walker takes the
realizer as an explicit argument. The string indirection carried a
real hazard class for zero payoff: a registry populated by import
side-effects is EMPTY in a fresh process until the right module
happens to be imported, a timing miss invisible to in-suite tests
(process-global state masks it). Dissolving the registry deleted
the hazard class instead of gating it.

The two extension axes the dual-registry design named survive
without it:

* Adding a new BC type means adding one
  :class:`BoundaryTraceLaw` subclass with a ``key=`` registration
  (+ the admission-table entries of the methods that support it).
  Each existing realizer adds a dispatch branch; no existing law
  changes.
* Adding a new transport method means minting one method-mesh
  (structurally conforming to ``TransportMethod``), one
  :class:`BoundaryRealizer`, and one admission table. **M** realizer
  branches need to be implemented (one per admitted BC), but no
  existing law changes — and no central registration step exists.


.. _bc-cross-method-stubs:

Adopting the architecture in a new method (and the retired stubs)
-----------------------------------------------------------------

A method adopts the unified BC architecture by shipping three
pieces (the diffusion adoption at #290 P3–P7b is the worked
example):

1. a **method-mesh** (``DiffusionMesh(MaterialMesh)``) carrying the
   method's trace + an admission table + the
   ``realize_boundary_law`` arm — it conforms structurally to
   :class:`~orpheus.transport.method.TransportMethod` and gets the
   whole tag → law → realized-``bc``-dict pipeline from the shared
   :func:`~orpheus.transport.method.resolve_boundary_conditions`
   body;
2. a **realizer** (``DiffusionBoundaryRealizer``) mapping each
   admitted law onto the method's operators;
3. a **method space** (``DiffusionMethodSpace``) carrying whatever
   discretization metadata the realizer reads.

.. note::

   **Historical: the Wave-5 stub scaffolding (retired #290 P7b).**
   Between Wave 5 and #290, ``MoCBoundaryRealizer`` /
   ``MCBoundaryRealizer`` / ``CPBoundaryRealizer`` (and, until #290
   P3, ``DiffusionBoundaryRealizer``) existed as
   ``NotImplementedError`` stubs auto-registered by each method's
   ``__init__.py``, "holding the dispatch architecture in place."
   With the registry dissolved there is no dispatch table to hold a
   place in — a stub realizer that cannot realize anything serves no
   consumer — so the three stub modules, their auto-import lines,
   and their stub-invariant tests were deleted. MoC / MC / CP keep
   their legacy solver-side BC validation until each adopts the
   architecture per the recipe above.

The SN realizer is **not** auto-imported by ``orpheus.sn.__init__``
(it's a heavy module that every SN consumer pays for); the
:class:`~orpheus.sn.mesh.augmented_mesh.SNMesh` construction imports
it explicitly when it needs it.


.. _bc-worked-example:

Worked example — end to end
===========================

The following walks the
``BC("vacuum") → VacuumInflow → SNBoundaryRealizer.realize →
IncomingOrdinateMaskTensor`` chain that
:meth:`orpheus.sn.mesh.augmented_mesh.SNMesh.realize_boundary_law`
performs per face (driven by the shared
:func:`~orpheus.transport.method.resolve_boundary_conditions` body) at
SNMesh construction time. The example uses a 1-D Cartesian slab; the same
chain runs on Mesh2D with face labels ``xmin`` / ``xmax`` /
``ymin`` / ``ymax``.

Step 1 — declaration on the mesh
--------------------------------

The user declares the vacuum BC on the mesh's left face:

.. code-block:: python

   from orpheus.geometry.mesh import Mesh1D, BC

   mesh = Mesh1D(
       edges=np.linspace(0.0, 1.0, 11),
       mat_ids=np.zeros(10, dtype=int),
       coord=CoordSystem.CARTESIAN,
       bc_left=BC("vacuum"),
       bc_right=BC("reflective"),
   )

The :class:`~orpheus.geometry.mesh.BC` dataclass is a thin wrapper
``BC(kind: str, params: dict)`` with no SN-specific knowledge. The
mesh is method-agnostic.

Step 2 — law resolution (in ``SNMesh.__init__``)
------------------------------------------------

When :class:`~orpheus.sn.mesh.augmented_mesh.SNMesh` is constructed against the
mesh, the shared
:func:`~orpheus.transport.method.resolve_boundary_conditions` body
walks the four (1-D: two) faces, calling
:meth:`~orpheus.sn.mesh.augmented_mesh.SNMesh.realize_boundary_law`
per face:

.. code-block:: python

   law_cls = SNMesh.BOUNDARY_OPERATOR_REGISTRY["vacuum"]
   # law_cls is VacuumInflow  (registry key -> law class lookup)
   law = law_cls()
   # law is a zero-arg instance: VacuumInflow has no parameters

The :attr:`SNMesh.BOUNDARY_OPERATOR_REGISTRY` is the SN-side view
of the law registry; today it carries only ``"vacuum"`` and
``"reflective"`` because those are the only kinds the SN sweep
pipeline has been wired for in production (the other three —
white, periodic, albedo — are realizable but require sweep-side
plumbing tracked in separate issues). Adding a new kind is one
dict-entry edit.

Step 3 — method space construction
----------------------------------

The per-face
:meth:`~orpheus.sn.mesh.augmented_mesh.SNMesh.realize_boundary_law`
calls share **one** unified
:class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace` for the whole
mesh, built once and stored on ``self._trace``. The factory is the
geometry-blind :meth:`AngularTraceSpace.from_quadrature_and_layout
<orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.from_quadrature_and_layout>`
— it takes the angular quadrature and the mesh's
:attr:`~orpheus.sn.mesh.augmented_mesh.SNMesh.boundary_face_layout` (a
:class:`~orpheus.numerics.face_layout.FaceLayout`, the single source of
truth for which faces exist and how they pack into one flat buffer), and
nothing else:

.. code-block:: python

   from orpheus.numerics.spaces.angular_trace_space import AngularTraceSpace

   self._trace = AngularTraceSpace.from_quadrature_and_layout(
       self.quad, self.boundary_face_layout,
   )

The trace stores the **signed projection** :math:`\Omega \cdot \hat
n_f` once per face as a ``(n_faces, N)`` float array — *not* two
direction-specific boolean masks. Inflow and outflow are
**selectors** over this one table, derived on demand by the sign of
:math:`\Omega \cdot \hat n_f` (the one-space-two-selectors design of
Issues #205 / #201; see :ref:`bc-trace-structure`). For the ``xmin``
face (``axis=0``, outward normal :math:`-\hat x`, so
:math:`\Omega \cdot \hat n = -\mu_x`), the inflow predicate
:math:`\Omega \cdot \hat n < -\epsilon` becomes :math:`-\mu_x[n] <
-\epsilon`, i.e. :math:`\mu_x[n] > \epsilon` — the rightward-pointing
ordinates are inflow at the left boundary, as expected.

The :meth:`SNMethodSpace.for_face` factory takes the precomputed trace
through a **single** ``trace=`` argument and extracts the per-face
inflow indices for the requested face:

.. code-block:: python

   from orpheus.sn.mesh.method_space import SNMethodSpace

   method_space = SNMethodSpace.for_face(
       mesh=self.mesh,           # optional metadata; None at d≥3
       quadrature=self.quad,
       face="xmin",
       trace=self._trace,
   )
   # for_face derives inflow_indices from the trace for this one face:
   #   inflow_indices = trace.inflow_indices_for_face("xmin")
   # i.e. the 1-D int array [n for n in range(N) if -mu_x[n] < -eps].

There is **one** trace object and **one** ``trace=`` parameter, not an
``inflow_trace`` / ``outflow_trace`` pair. The directional split lives
in the two selector methods
(:meth:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.inflow_indices_for_face`
/
:meth:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.outflow_indices_for_face`),
not in the space's identity. Why one space? Because whether an ordinate
is incoming or outgoing at a face is a *predicate* evaluated against the
same boundary data, not a property of two distinct domains — folding the
two former spaces into one removes a class of bugs where the inflow and
outflow descriptions of the *same* boundary could drift out of sync, and
gives the Wave-O adjoint work (#208) a single
:math:`|\Omega\cdot\hat n|`-weighted boundary inner product to install
(see :ref:`bc-trace-structure`). The ``mesh=`` argument is optional
metadata: nothing in the realizer chain reads it (the inflow indices
come from the trace), so an axis-native ``SNMesh`` with no legacy mesh
adapter passes ``None`` (C5.3, #225).

.. note:: **Historical — the pre-Issue-#188 split trace.** Before #188
   and the #205 / #201 unification, this step built a per-face
   ``InflowTraceSpace`` / ``OutflowTraceSpace`` *pair* via
   ``InflowTraceSpace.from_mesh_and_quadrature(mesh, quad,
   faces=("left", "right"))``, and :meth:`SNMethodSpace.for_face` took
   *two* arguments, ``inflow_trace=`` and ``outflow_trace=``. That
   machinery is retired: there is now one geometry-blind
   :class:`AngularTraceSpace <orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace>`
   with two directional selectors. This note records the older API only
   so that references to ``_inflow_trace`` / ``_outflow_trace`` in
   pre-#188 commit history are legible; it is **not** the current code.

Step 4 — realization
--------------------

The :class:`SNBoundaryRealizer` is now invoked — directly, by the
method-mesh that owns it (``SNMesh.realize_boundary_law``, the SN
arm of the :class:`~orpheus.transport.method.TransportMethod` hook;
#290 P7b removed the registry lookup that used to sit here).
Instantiation is stateless:

.. code-block:: python

   from orpheus.sn.boundary.realizer import SNBoundaryRealizer

   realized = SNBoundaryRealizer().realize(law, method_space)

The realizer's vacuum branch fires:

.. code-block:: python

   # Inside SNBoundaryRealizer.realize:
   if isinstance(law, VacuumInflow):
       return IncomingOrdinateMaskTensor(
           inflow_indices=method_space.inflow_indices,
           n_ordinates=quad.N,
           axis=0,
       )

The returned ``realized`` is a
:class:`~orpheus.numerics.operator.IncomingOrdinateMaskTensor`,
which is a self-adjoint, idempotent Wave-0 primitive that reports
``is_adjointable = True`` (:math:`M = M^{\mathsf T}`) and
``is_invertible = False`` (a rank-deficient projection — no inverse).

Step 5 — wrap with a kind tag
-----------------------------

Every ``SNMesh.bc[<face>]`` entry carries a uniform 1-arg
``apply(psi)`` contract (Wave 9 migrated 13 production sites from
2-arg to 1-arg; C4 / #220 re-keyed the per-attribute ``bc_<face>``
surface to the face-name-keyed :attr:`~SNMesh.bc` dict — see
:ref:`bc-face-name-carve`). Post Issue #186 / C-B3.4 the
:class:`~orpheus.geometry.boundary._bound_compat._BoundBoundaryOperator`
shim is a **strict 1-arg passthrough** that adds two pieces of
metadata to the realized operator:

* a free-form ``kind`` string tag — since #290 P7b read off the
  realized LAW's own registry key (``law.key``; identical to the
  declared :class:`~orpheus.geometry.mesh.BC` kind because every
  admission-table entry maps a tag to the law registered under that
  same key) — load-bearing for the
  ``sn_mesh.bc["xmin"] == "vacuum"`` string-equality surface that
  several SN tests rely on;
* :attr:`~orpheus.numerics.operator.LinearOperator.is_invertible` /
  :attr:`~orpheus.numerics.operator.LinearOperator.is_adjointable`
  delegation to the wrapped inner operator so consumers composing the
  shim with other Wave-0 primitives inherit the right surface.

The shim's ``apply`` / ``apply_transpose`` signatures are strict
1-arg ``(self, psi)`` — extra positional or keyword arguments
raise :class:`TypeError`. The pre-Issue-#186 affordance that
swallowed ``*_extra, **_kw`` was the last remnant of the 2-arg
legacy era; it was dropped in C-B3.4 alongside the descriptor
cleanup because every production and test call site is now
strict 1-arg.

.. code-block:: python

   from orpheus.geometry.boundary._bound_compat import _BoundBoundaryOperator

   return _BoundBoundaryOperator(realized, kind=law.key)

The shim is **internal** to the package (not in :attr:`__all__`)
— a test pins its private status.

**Historical note (pre Issue #176 / #186).** The Wave-8/9
implementation carried an optional ``quadrature=`` kwarg that,
when non-``None``, bound an
``AngularQuadrature`` and forwarded
``inner.apply(psi, bound_quad)`` to a legacy 2-arg
:class:`BoundaryTraceLaw` body. That dual-mode existed ONLY
because the trace factory (then named ``InflowTraceSpace.from_mesh_and_quadrature``,
since C5.3 the geometry-blind
:meth:`AngularTraceSpace.from_quadrature_and_layout
<orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.from_quadrature_and_layout>`)
raised :class:`NotImplementedError` for curvilinear ``Mesh1D``, which
forced the per-face resolution (then ``SNMesh._resolve_one``, since
#290 P7b ``SNMesh.realize_boundary_law``) to bypass the realizer for
spherical / cylindrical meshes. Issue #188 lifted that
deferral; Issue #176 (C176.1) dropped the bound-quadrature mode
here because no production-issued shim carried
``_quadrature is not None`` after C188.3. Issue #186 (C-B3.4)
then dropped the residual ``*_extra, **_kw`` argument-swallow
because, with concrete-BC ``apply`` methods retired and all
production / test sites strict 1-arg, the defensive net was
dead code.

Step 6 — consumption by the sweep
---------------------------------

At each sweep call site the resolved operator is applied with
the uniform 1-arg interface:

.. code-block:: python

   psi_in = self.bc["xmin"].apply(psi_out_full)
   # Shim forwards to IncomingOrdinateMaskTensor.apply, which zeros
   # only the inflow ordinates; outflow rows pass through unchanged.

The downstream sweep reads ``psi_in[inflow_ord]`` only — every
production call site was audited in Wave 8 (see the §16A.5
compatibility-audit table in the Wave 8 close-out memo and the
:ref:`bc-vacuum-semantic-correction` section below).


.. _bc-universal-invariants:

Universal ``assert_*`` invariants
=================================

The :class:`~orpheus.geometry.boundary.BoundaryTraceLaw` ABC ships
**five universal** assertion methods (no-op defaults; concrete laws
override) plus **four specific** assertions on the BCs that need
them. Together they form the structural verification surface that
the :mod:`tests.geometry.test_bc_universal_invariants` suite
exercises (~30 L1 tests).

The split between "universal" and "specific" follows the
Grand Report v3 §16A.12 + §27.6 catalog. Universal invariants are
properties the affine law :eq:`affine-bc-form` should satisfy for
**any** physically meaningful BC; specific invariants pin properties
that only a subset of laws claim (e.g. involution is meaningful for
reflective but not for white).

.. list-table:: Universal invariants
   :header-rows: 1
   :widths: 30 25 45

   * - Method
     - Pinned error
     - What it asserts
   * - ``assert_inflow_outflow_classification``
     - :class:`~orpheus.geometry.boundary.IncomingOutgoingTraceClassificationError`
       (ERR-040)
     - Every ordinate at the face is either inflow or outflow (no
       tangential ordinates allowed by the law's contract). Default:
       no-op; reflective overrides to require strict partition.
   * - ``assert_outgoing_leakage_unconstrained``
     - n/a (architectural contract)
     - The outgoing trace flux is not constrained by the BC.
       Default: no-op. Future Dirichlet-outflow / prescribed
       cell-edge interface laws would override.
   * - ``assert_geometry_map_measure_preserving``
     - :class:`~orpheus.geometry.boundary.BoundaryGeometryMapNotMeasurePreservingError`
       (ERR-042)
     - The geometric map :math:`G` preserves the angular measure
       :math:`w(\Omega)\,|\Omega \cdot \hat n|`. Default: no-op.
       Reflective overrides (delegates to the involution check).
   * - ``assert_response_positive_if_declared``
     - :class:`~orpheus.geometry.boundary.BoundaryResponseNotPositiveError`
       (ERR-043)
     - If a response kernel is declared, it produces non-negative
       output on the inflow trace. Default: no-op. Albedo / white
       override.
   * - ``assert_source_lives_on_incoming_trace``
     - :class:`~orpheus.geometry.boundary.BoundarySourceNotOnIncomingTraceError`
       (ERR-047)
     - The source :math:`q` is nonzero only on :math:`\Gamma_-`.
       Default: no-op (the :class:`NoSource` default trivially
       satisfies). Prescribed-inflow + future user-source classes
       override.

.. list-table:: Specific invariants
   :header-rows: 1
   :widths: 40 25 35

   * - Method
     - Pinned error
     - Where it's defined
   * - :meth:`ReflectiveBoundary.assert_is_involutive`
     - :class:`~orpheus.geometry.boundary.ReflectionNotInvolutiveError`
       (ERR-044)
     - :class:`ReflectiveBoundary` (reflection-index table must
       satisfy ``perm[perm] == arange``).
   * - :meth:`ReflectiveBoundary.assert_maps_inflow_to_outflow`
     - :class:`~orpheus.geometry.boundary.ReflectionDidNotMapInflowToOutflowError`
       (ERR-045)
     - :class:`ReflectiveBoundary` (every inflow ordinate maps to an
       outflow ordinate under the reflection).
   * - :meth:`ReflectiveBoundary.assert_direction_norm_preserved`
     - inherits ERR-042 via measure-preservation
     - :class:`ReflectiveBoundary` (delegates to involution check).
   * - :meth:`WhiteBoundary.assert_submarkov` /
       :meth:`AlbedoBoundary.assert_submarkov`
     - :class:`~orpheus.geometry.boundary.SubmarkovViolationError`
       (ERR-046)
     - :class:`WhiteBoundary` + :class:`AlbedoBoundary` (the sub-Markov
       kernel constraint :math:`\int R\,\mathrm{d}y \leq 1` per row;
       albedo :math:`> 1` violates this physically).


.. _bc-named-error-catalog:

Named-error catalog (ERR-040..ERR-047)
======================================

Per Grand Report v3 §26A.4 and the `vv-principles` skill's
"Log every caught bug" directive, every typed error has a matching
``@pytest.mark.catches("ERR-NNN")`` decorator on the test that
proves it fires under the right fault-injection. The eight errors
shipped under :mod:`orpheus.geometry.boundary._errors` are:

.. list-table::
   :header-rows: 1
   :widths: 18 26 12 44

   * - Error class
     - Trigger condition
     - Mode
     - Mechanism
   * - :class:`~orpheus.geometry.boundary.IncomingOutgoingTraceClassificationError`
       (ERR-040)
     - Tangential ordinate at a face that requires strict partition.
     - #5 (index)
     - ``assert_inflow_outflow_classification`` finds
       ``|Ω · n| ≤ ε`` on a face where the law's contract
       forbids it.
   * - :class:`~orpheus.geometry.boundary.VacuumAppliedToOutgoingTraceError`
       (ERR-041)
     - Vacuum law applied to an outgoing trace.
     - #6 (convention)
     - Vacuum sets only :math:`\gamma_- \psi = 0`; applying it on
       :math:`\Gamma_+` is geometrically meaningless and typically
       indicates a wrong face annotation.
   * - :class:`~orpheus.geometry.boundary.BoundaryGeometryMapNotMeasurePreservingError`
       (ERR-042)
     - Geometric map :math:`G` does not preserve
       :math:`w(\Omega) |\Omega \cdot \hat n|`.
     - #5 + #6
     - Wrong reflection-index table or inconsistent quadrature
       :math:`\mu_n` / weights.
   * - :class:`~orpheus.geometry.boundary.BoundaryResponseNotPositiveError`
       (ERR-043)
     - Response kernel produces negative output.
     - #1 (sign)
     - Sign-flipped kernel construction (e.g. accidental
       :math:`-\alpha`).
   * - :class:`~orpheus.geometry.boundary.ReflectionNotInvolutiveError`
       (ERR-044)
     - Reflection permutation is not its own inverse:
       :math:`\pi \circ \pi \neq \mathrm{id}`.
     - #5 (index)
     - Wrong reflection axis or partial permutation in the
       :meth:`quadrature.reflection_index` table.
   * - :class:`~orpheus.geometry.boundary.ReflectionDidNotMapInflowToOutflowError`
       (ERR-045)
     - Reflection maps an inflow ordinate to itself.
     - #5 (index)
     - Non-axis-aligned reflection mislabeled as ``ReflectiveBoundary``
       (the right BC family would be a future ``SymmetryBoundary``).
   * - :class:`~orpheus.geometry.boundary.SubmarkovViolationError`
       (ERR-046)
     - Sub-Markov BC with :math:`\alpha > 1`.
     - #4 (factor)
     - Albedo / white kernel scalar exceeds 1.0 — physically this
       would imply a source on the boundary surface.
   * - :class:`~orpheus.geometry.boundary.BoundarySourceNotOnIncomingTraceError`
       (ERR-047)
     - Boundary source :math:`q` has nonzero outflow entries.
     - #6 (convention)
     - User-supplied :class:`InflowSourceSpec` has nonzero entries on
       :math:`\Gamma_+`; geometrically meaningless and indicates a
       wrong source-shape contract.

All eight extend :class:`ValueError` (via the
:class:`~orpheus.geometry.boundary.BoundaryError` base) so existing
``except ValueError`` consumers from the pre-refactor code keep
working. Every catch site can additionally pattern-match on the
typed subclass to recover the offending law name from
:attr:`BoundaryError.law`.


.. _bc-rank-n-algebra:

.. _bc-descriptor-tree-vs-operator-tree:

Descriptor-tree algebra for rank-N boundaries
=============================================

Rank-N (Marshak, partial-current) boundary conditions are
**not** a special class. They are expressed via a closed
**descriptor-tree algebra** over
``BoundaryTraceLaw | LawSum | LawScaled`` nodes — pure
declarative structure with **no** ``apply`` method on any node.
The :class:`~orpheus.geometry.boundary.BoundaryTraceLaw` algebra
dunders (``+``, ``-``, ``*``, ``/``, unary ``-``) return
:class:`~orpheus.geometry.boundary.LawSum` /
:class:`~orpheus.geometry.boundary.LawScaled` instances, never
operators. The :func:`~orpheus.geometry.boundary.realize_recursively`
type transformer (method-blind since #290 P7b — the leaf realizer is
an explicit argument) is the **sole** path from descriptor tree to
operator tree.

The §15.2 sum-of-tensor-products form

.. math::
   :label: bc-rank-n-tensor-decomposition

   R \;=\; \sum_{\alpha} c_{\alpha}\, G_{\alpha},
   \qquad c_{\alpha} \in \mathbb{R},
   \quad G_{\alpha} \in
   \{\text{permutation, average, mask, wrap, identity, source}\},

.. vv-status: bc-rank-n-tensor-decomposition documented

maps onto the LawXxx algebra as ``c1 * law_1 + c2 * law_2 + ...``,
where each ``c_i * law_i`` term is a :class:`LawScaled` node and
the sum is a :class:`LawSum` node.

The standard Marshak boundary (Bell & Glasstone 1970 §1.5) — a
mix of specular reflection (weight :math:`c_1`) and diffuse
white reflection (weight :math:`c_2`) — is:

.. code-block:: python

   from orpheus.geometry.boundary import (
       LawScaled, LawSum,
       ReflectiveBoundary, WhiteBoundary,
       realize_recursively,
   )
   from orpheus.sn.boundary.realizer import SNBoundaryRealizer
   from orpheus.sn.mesh.method_space import SNMethodSpace

   # Build the descriptor tree — no realization yet.
   spec = ReflectiveBoundary(axis="x", albedo=1.0)
   white = WhiteBoundary(axis="x", outward_sign=+1, albedo=1.0)
   marshak_law = 0.3 * spec + 0.7 * white
   # marshak_law is:
   #   LawSum(
   #       LawScaled(0.3, ReflectiveBoundary(axis="x", albedo=1.0)),
   #       LawScaled(0.7, WhiteBoundary(axis="x", outward_sign=+1,
   #                                    albedo=1.0)),
   #   )
   # NOT callable — no .apply method on LawSum or its children.
   assert not hasattr(marshak_law, "apply")

   # Realize the tree at one face. realize_recursively walks
   # LawSum / LawScaled / leaf-law nodes and emits the matching
   # Wave-0 operator-tree composers around realized 1-arg leaves.
   # The walker is method-blind — pass the method's own realizer.
   ms = SNMethodSpace.for_face(...)
   marshak_op = realize_recursively(marshak_law, ms, SNBoundaryRealizer())
   # marshak_op is:
   #   OperatorSum(
   #       ScaledOperator(0.3, PermutationOperator(...)),
   #       ScaledOperator(0.7, AngularAverageOperator(...))
   #   )
   psi_in = marshak_op.apply(psi_out)   # 1-arg LinearOperator

The output is a Wave-0
:class:`~orpheus.numerics.operator.OperatorSum` of
:class:`~orpheus.numerics.operator.ScaledOperator`-wrapped Wave-0
primitives, consumable by the SN sweep / Krylov path via the
uniform 1-arg :meth:`apply`. The descriptor-tree algebra and the
operator-tree algebra are **separate type families**: the
descriptor tree is built with ``LawXxx`` nodes that have **no**
``apply``; the operator tree is built with ``OperatorXxx`` nodes
that **do** have ``apply``. The two families never inter-compose —
mixing a :class:`LawNode` with an already-realized
:class:`~orpheus.numerics.operator.LinearOperator` via ``+`` is a
type error; the user MUST call :func:`realize_recursively` first.

Closed-algebra guarantees
-------------------------

* **Constant folding on scalars.**
  ``LawScaled(α, LawScaled(β, x))`` collapses to
  ``LawScaled(α * β, x)`` at construction time. The intermediate
  ``LawScaled`` nesting never appears at rest, which keeps the tree
  shallow under repeated scalar multiplication.
* **No associativity flattening on sums.** ``(a + b) + c`` is
  :class:`LawSum(LawSum(a, b), c)`, distinct from
  :class:`LawSum(a, LawSum(b, c))`. The walker treats both shapes
  identically — the realized output is the same Wave-0 operator
  algebra value up to floating-point non-associativity in the
  final sum.
* **Subtraction rewrites via :class:`LawScaled(-1, ...)`.** The
  unary ``-`` operator and the binary ``-`` operator both produce
  trees containing only :class:`LawSum` / :class:`LawScaled`
  nodes — there is no dedicated ``LawDifference`` type.
* **Division rewrites via :class:`LawScaled(1/α, ...)`.** Pure
  syntactic sugar for ``LawScaled(α, ...).__truediv__``.

The pre-refactor implementations
--------------------------------

Two prior approaches converged on the present descriptor-tree
design through empirical falsification:

**Wave 11 (~2026-03)** — ``MixedBoundaryOperator(components:
list[tuple[float, BoundaryOperator]])`` class whose
:meth:`apply` body looped over ``components`` and summed
``coeff * primitive.apply(psi, quad)``. The SN realizer
dispatched on it via an ``isinstance(law, MixedBoundaryOperator)``
branch that ran the same loop with
``coeff * realize(primitive, ms)`` summed via
:class:`OperatorSum`. Wave 11 deleted this code because the
delayed-realization-by-container pattern broke down once vacuum
needed per-face inflow indices that the bare-law container had
no access to.

**β1 interim landing (Issue #186 / B3, ~2026-04)** — every
:class:`BoundaryTraceLaw` inherited the Wave-0 operator-algebra
dunders from :class:`~orpheus.numerics.operator.LinearOperator`,
so writing ``0.3 * spec + 0.7 * white`` directly produced an
:class:`OperatorSum` of :class:`ScaledOperator`-wrapped raw
:class:`BoundaryTraceLaw` leaves (NOT realized). The
:func:`realize_recursively` walker then traversed the Wave-0
composer tree, realized each leaf, and emitted a parallel
operator tree. This achieved the right algebraic shape but
**violated the type system**: the resulting tree was an
:class:`OperatorSum` instance (a :class:`LinearOperator`!) whose
:meth:`apply` could not actually be called before realization —
calling it raised :class:`BoundaryError` at apply-time because the
leaves were laws, not operators. The convention "you must realize
this OperatorSum before calling apply" was a runtime contract that
the type system did nothing to enforce. β1 retained the
ergonomic of "the same ``+`` operator before and after
realization" but at the cost of conflating two type families.

**β2 (this scope, Issue #186 / B3 + β2)** — separates the two
type families explicitly: :class:`LawSum` / :class:`LawScaled` for
the descriptor tree (no :meth:`apply`); :class:`OperatorSum` /
:class:`ScaledOperator` for the operator tree (with
:meth:`apply`). The static type system enforces "you cannot call
this until it's realized" — :class:`LawSum` has no :meth:`apply`
method on the class, so the linter flags ``tree.apply(...)``
without running the program. The ergonomic of "the same ``+``"
survives because both the law-tree and the operator-tree use the
same Python ``+`` syntax; the runtime dispatch on type tells the
reader (and the type checker) which algebra is in effect.

The Wave 6 snapshot harness verifies that the bit-identity of the
Marshak case ``0.3 * spec + 0.7 * white`` is **preserved** by the
β2 transition: the realized
``OperatorSum(ScaledOperator, ScaledOperator)`` reduction tree
matches the β1-era output to the same ULP tolerance, because the
operator-tree shape after :func:`realize_recursively` is
algebraically identical.


.. _bc-realize-recursively:

The ``realize_recursively`` walker — a descriptor → operator type transformer
=============================================================================

:func:`~orpheus.geometry.boundary.realize_recursively` is the
**type transformer** from the descriptor-tree algebra
(``BoundaryTraceLaw | LawSum | LawScaled``) to the operator-tree
algebra (``LinearOperator`` with
:class:`~orpheus.numerics.operator.OperatorSum` /
:class:`~orpheus.numerics.operator.ScaledOperator` composers
around realized 1-arg leaves). Calling it is the **only** path
from a non-callable descriptor to a callable operator.

The dispatch is exhaustive on the descriptor-tree node types:

.. code-block:: python

   def realize_recursively(
       node: BoundaryTraceLaw | LawSum | LawScaled,
       method_space: SNMethodSpace,
   ) -> LinearOperator:
       if isinstance(node, BoundaryTraceLaw):
           # Leaf: realize via the SN realizer registry.
           return SNBoundaryRealizer().realize(node, method_space)
       if isinstance(node, LawScaled):
           # Scalar-times-law: wrap the realized inner in ScaledOperator.
           inner_op = realize_recursively(node.inner, method_space)
           return ScaledOperator(node.scalar, inner_op)
       if isinstance(node, LawSum):
           # Sum: realize each side, wrap in OperatorSum.
           a_op = realize_recursively(node.a, method_space)
           b_op = realize_recursively(node.b, method_space)
           return OperatorSum(a_op, b_op)
       raise TypeError(
           f"realize_recursively expected BoundaryTraceLaw | LawSum | "
           f"LawScaled, got {type(node).__name__}."
       )

Usage on the descriptor tree:

.. code-block:: python

   from orpheus.geometry.boundary import (
       ReflectiveBoundary, WhiteBoundary, realize_recursively,
   )
   from orpheus.sn.boundary.realizer import SNBoundaryRealizer

   # Build the descriptor tree.
   law = (
       0.3 * ReflectiveBoundary(axis="x")
       + 0.7 * WhiteBoundary(axis="x", outward_sign=+1)
   )
   # law is LawSum(LawScaled(0.3, ...), LawScaled(0.7, ...)).

   # Realize once, at face resolution time — the walker is
   # method-blind, so the method's realizer is passed explicitly.
   realized = realize_recursively(law, method_space, SNBoundaryRealizer())
   # realized is:
   #   OperatorSum(
   #       ScaledOperator(0.3, PermutationOperator(...)),
   #       ScaledOperator(0.7, AngularAverageOperator(...)),
   #   )
   psi_in = realized.apply(psi_out)   # 1-arg LinearOperator

Type-system contract
--------------------

The walker's input is intentionally narrow:

* :class:`BoundaryTraceLaw` instances and the two descriptor-tree
  composer dataclasses (:class:`LawSum`, :class:`LawScaled`) are
  the only valid node shapes.
* Wave-0 operator-tree composers
  (:class:`~orpheus.numerics.operator.OperatorProduct`,
  :class:`~orpheus.numerics.operator.TensorProductOperator`,
  :class:`~orpheus.numerics.operator.SumOfTensorProductsOperator`,
  :class:`~orpheus.numerics.operator.OperatorSum`,
  :class:`~orpheus.numerics.operator.ScaledOperator`) are **not**
  recognized — they belong to the operator tree, not the
  descriptor tree, so they should never appear in the realizer's
  input.
* Unknown nodes raise :class:`TypeError` (not
  :class:`BoundaryError`) with the offending type name in the
  message, because this is a **typing** failure (caller passed
  the wrong kind of object), not a BC-domain failure.

The β1 → β2 transition (see :ref:`bc-rank-n-algebra`) sharpened
the walker's type signature from "any Wave-0 composer tree with
:class:`BoundaryTraceLaw` leaves" to "any descriptor-tree node".
The dispatch table shrank from five Wave-0 composers + leaf to
three descriptor types + leaf (counting the leaf as the same
:class:`BoundaryTraceLaw` branch). The eliminated branches
(:class:`OperatorProduct`, :class:`TensorProductOperator`,
:class:`SumOfTensorProductsOperator`) handle operator composition
patterns that have no descriptor-tree analog — they were dead
dispatch paths once :class:`LawSum` / :class:`LawScaled` replaced
the in-tree Wave-0 algebra. Removing them clarified the walker's
role: it is **exactly** the type transformer between the two
algebras, nothing more.

Placement — the deferral that fired at the second method
--------------------------------------------------------

The walker lives in :mod:`orpheus.geometry.boundary` (the
``_realizer`` module — the realization seam, next to the
:class:`~orpheus.geometry.boundary.BoundaryRealizer` Protocol it
dispatches through), and it is **method-blind**: the leaf realizer
is a REQUIRED argument, and ``method_space`` is threaded verbatim to
that realizer without inspection. It is *not* on the production
single-BC path — production realizes one BC directly (each
method-mesh's ``realize_boundary_law`` arm →
``<Method>BoundaryRealizer().realize``). The walker is the **rank-N
composition entry point**: the only thing that realizes a
*descriptor tree* (the Marshak ``0.3 * Reflective + 0.7 * White``
partial-current BC of :eq:`affine-bc-form`) rather than a single
leaf law. Production does not yet wire a rank-N BC, so the walker
runs only from the rank-N tests.

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
on the **method-mesh layer** (``SNMesh`` / ``DiffusionMesh`` — the
user ruling: the method-mesh IS the method's behavior carrier, not a
stateless singleton), the twin per-mesh ``_resolve_bcs`` loops
collapsed into the ONE shared
:func:`~orpheus.transport.method.resolve_boundary_conditions` body,
and the walker moved here with an explicit-realizer signature. Only
the **leaf** dispatch was ever method-bound; the **composer**
dispatch (:class:`LawSum` / :class:`LawScaled` → operator composers)
survived the move byte-identically, exactly as recorded.


.. _bc-vacuum-semantic-correction:

The vacuum semantic correction (§16A.5)
=======================================

The most subtle design decision of the refactor concerns vacuum.
Pre-Wave-7 the legacy ``VacuumBoundaryOperator.apply`` body was:

.. code-block:: python

   def apply(self, psi_out: np.ndarray, quadrature) -> np.ndarray:
       return np.zeros_like(psi_out)

This returns **all zeros**, including the outflow ordinates that
the BC has no physical interpretation for (vacuum says nothing
about :math:`\gamma_+ \psi`; it only sets :math:`\gamma_- \psi = 0`).

The post-Wave-8 SN realizer's vacuum branch returns
:class:`~orpheus.numerics.operator.IncomingOrdinateMaskTensor` —
a sparse mask that zeros **only the inflow ordinates** and
preserves the outflow trace. This is the §16A.10 trace-correct
representation: the operator's apply is a projector onto the
inflow ordinate subspace, which is the right algebraic object
for the affine law :eq:`affine-bc-form` to read.

Why this matters
----------------

Three downstream consequences make the inflow-only mask the right
contract:

1. **Sensitivity adjoints.** A future adjoint sensitivity path
   needs the outflow trace preserved to compute the response
   of an outgoing-current functional to the inflow BC. The
   zeros-all body would silently lose the gradient at the
   boundary.
2. **Compositional clarity.** The realized vacuum operator is
   now self-adjoint and idempotent (it's a projector). The
   zeros-all operator is also idempotent but is the
   ``ZeroOperator`` projector, which is the wrong type tag
   for "inflow-mask only".
3. **Algebraic uniformity.** Every other rank-1 law (reflective,
   white, albedo) acts on the **inflow ordinates only** and
   leaves the outflow rows untouched. The legacy vacuum was
   the asymmetric special case; the post-Wave-8 vacuum joins
   the family.

The algebra: where the two semantics agree and where they diverge
-------------------------------------------------------------------

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

.. (vv-status rationale) Structural/explanatory identity: contrasts the legacy
   zeros-all vacuum with the trace-correct inflow-only mask (they agree on the
   inflow set, diverge on the outflow set). The trace-correct inflow-only
   semantics are pinned by the realizer snapshot gates
   (``tests/geometry/test_bc_equivalence_snapshot.py``, the §16A.5 semantic
   capture). An explanatory comparison, not a separate solver claim.
.. vv-status: vacuum-legacy-vs-trace-correct documented

The two functions agree **on** :math:`I_f` (both give 0) and
**diverge** on :math:`O_f` (legacy gives 0, trace-correct gives
:math:`\psi_{\text{out}}[n]`). They diverge on :math:`T_f` too
in principle, but ORPHEUS's quadrature adapters carry every
tangential ordinate at :math:`\mu = 0` so :math:`\psi_{\text{out}}[n] = 0`
on :math:`T_f` for a properly-initialised flux — making the
divergence physically restricted to :math:`O_f`.

The §16A.5 production-relevant subset is **the inflow rows**.
Every SN sweep call site reads :math:`\psi_{\text{in}}[n]` only for
:math:`n \in I_f` — outflow rows are never consumed downstream.
The Wave 8 close-out audited all 13 ``bc.apply(...)`` sites in
``orpheus.sn.loss_representation`` (the dissolved ``sweep.py``) and :mod:`orpheus.sn.operators.streaming`:

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

No call site reads :math:`\psi_{\text{in}}[n]` for
:math:`n \in O_f`, so the two semantics produce **bit-identical
observable output** for every existing production consumer.
This is what makes the §16A.5 trace-correct representation a
safe semantic upgrade: it strictly extends the correctness
boundary (outflow rows are now correct under the new algebra
where they were ill-defined under the old one) without changing
any observable downstream value.

Post Issue #188 + #176 (2026-05-11) the **realizer path is
uniform** across every supported mesh — 1-D Cartesian, 1-D
spherical, 1-D cylindrical, and 2-D Cartesian. The curvilinear
sweeps now consume the realizer-routed
:class:`IncomingOrdinateMaskTensor` exactly like the slab sweep
does. **Empirical confirmation**: spherical 26/26 + cylindrical
25/25 + MMS curvilinear 2/2 xfail (pre-existing ERR-026) green
on C188.3 — the curvilinear sweeps were previously consuming the
zeros-all legacy body via the bound-quadrature shim, and now
consume the inflow-only-mask realizer output. Production result
is bit-identical on inflow rows (the only rows the sweep reads),
confirming the Wave 8 call-site audit empirically.

The Wave 6 snapshot harness gates this explicitly: the
``vacuum_lebedev17`` case compares **inflow-row outputs only**,
documenting the intentional semantic divergence in a comment
adjacent to the test case.

"Option a" (Wave 7) — historical context, retired Issue #186
-------------------------------------------------------------

The Wave 7 brief considered three migration strategies for the
2-arg legacy path. **Option (a)** ("vacuum-stays-legacy") landed:
:class:`VacuumInflow` carried a standalone
:meth:`apply(psi_out, quad)` whose body was
``np.zeros_like(psi_out)`` (the pre-§16A.5 zeros-all form), and
the realizer path produced the inflow-only-mask form via
:class:`~orpheus.numerics.operator.IncomingOrdinateMaskTensor`.
The two paths agreed on inflow rows (the production-relevant
subset) and diverged on outflow rows.

Options (b) ("face-aware BC" — add ``face`` to the constructor)
and (c) ("combined ABC merge") were rejected because they would
have distorted the law's interface to better serve a transitional
path the refactor was retiring anyway.

**Status after Issue #186 (2026-05-11):** Option (a)'s standalone
``apply`` body has been **deleted**. :class:`VacuumInflow` (like
every other concrete BC) is a pure descriptor with no ``apply``
method; the only path to vacuum action is
:func:`realize_recursively` or :class:`SNBoundaryRealizer`
producing :class:`IncomingOrdinateMaskTensor` output. The
**two-paths-divergence is therefore eliminated by design** — there
is no longer a "second path" that could disagree with the realizer
path. The §16A.5 inflow-only-mask body is the **unique** vacuum
semantics in the post-#186 codebase.

This is the load-bearing architectural payoff of B3 + β2: the
documentation no longer needs to caveat which path you're on,
because there is only one path. See
:ref:`bc-trace-law-descriptor-model` for the design rationale.


.. _bc-trace-law-descriptor-model:

The trace-law descriptor model
==============================

Issue #186 / Scope B3 + β2 (landed 2026-05-11 on branch
``feature/bc-curvilinear-realizer-cleanup``) is the architectural
**closure** of the BC trace-law refactor. It collapses the
remaining 2-arg ``apply`` affordance from the Wave-8/9 era into a
**pure-descriptor** contract:

* :class:`BoundaryTraceLaw` no longer inherits
  :class:`~orpheus.numerics.operator.LinearOperator`.
* The abstract :meth:`apply` method that the mixin used to provide
  is gone. So is ``apply_transpose``. So is any operator-surface
  advertisement — the two-axis ``is_invertible`` / ``is_adjointable``
  predicates (before the #226 carve P4, the retired
  ``capabilities: ClassVar[frozenset[str]]`` frozenset).
* Every concrete BC (vacuum / reflective / white / albedo /
  periodic / prescribed-inflow) is now a **frozen dataclass**
  carrying only its parameters (axis, albedo, source, ...), its
  :attr:`kind` tag, its :attr:`creates_sweep_cycle` ``ClassVar``,
  its :attr:`geometry_map` / :attr:`response_kernel` /
  :attr:`source` property overrides, and the relevant
  :meth:`assert_*` invariants. **No** ``apply`` method on any
  concrete BC.
* The base class :class:`BoundaryTraceLaw` carries a **minimal
  algebra** that returns :class:`LawSum` / :class:`LawScaled`
  nodes — the descriptor-tree composition algebra documented at
  :ref:`bc-rank-n-algebra`. The dunders are: ``+``, ``-``, ``*``,
  ``/``, unary ``-``, plus their reflected variants. Each returns
  a new descriptor-tree node; none returns an operator.
* :class:`~orpheus.sn.boundary.realizer.SNBoundaryRealizer` is the
  **sole** bridge from descriptor to callable. There is no
  alternative path. Calling ``law.apply(psi)`` raises
  :class:`AttributeError` at runtime; a static type checker flags
  the call without running it.

The §16A.3 three-layer architecture (descriptor / realizer /
operator) is now **enforced by the type system**, not by convention.

What was tried and rejected before B3 + β2 landed
-------------------------------------------------

This subsection preserves the design history because the rejected
paths are the load-bearing intellectual content of the close-out
— future sessions asking "why does the realizer exist?" need to
see why every alternative failed.

**Option A** (Issue #176 / C176.3, ~2026-04). Concrete BC
:meth:`apply` methods kept a keyword-optional
``quadrature: AngularQuadrature | None = None`` parameter with
defensive errors:

* :class:`ReflectiveBoundary.apply` / :class:`WhiteBoundary.apply`
  raised :class:`BoundaryError` when ``quadrature is None``
  because their geometric / response operators needed the
  quadrature to construct themselves.
* :class:`VacuumInflow` / :class:`AlbedoBoundary` /
  :class:`PeriodicBoundary` / :class:`PrescribedInflow` accepted
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
   :class:`~orpheus.numerics.operator.IncomingOrdinateMaskTensor`
   output (the §16A.5 inflow-only mask). The two paths agreed on
   inflow rows (the production-relevant subset; see
   :ref:`bc-vacuum-semantic-correction`) but diverged on outflow
   rows. The divergence was harmless at every existing production
   call site but was a documentation-burden landmine for future
   adjoint-sensitivity consumers that read outflow rows.
3. **Liskov violation.** The abstract
   :meth:`BoundaryTraceLaw.apply(self, psi_out)` was strict 1-arg
   (post Issue #176 / C176.4). The concrete
   :meth:`apply(self, psi_out, quadrature=None)` was technically
   Liskov-substitutable (the optional parameter has a default), but
   calling ``bc.apply(psi)`` polymorphically on a
   :class:`BoundaryTraceLaw`-typed parameter could fail at runtime
   for Reflective/White — the type signature said "this works" but
   the runtime behaviour said "this raises". The static type
   system could not catch the failure because the contract was
   carried only in the docstring.

The B3 audit cataloged every remaining 2-arg ``bc.apply(psi,
quad)`` call site in production AND tests; none was load-bearing
for correctness. The Wave-6 snapshot regression carried legacy
halves (regenerated through the realizer path), the
:class:`PrescribedInflow` invariant tests ignored the
``quadrature`` argument anyway, and the realizer-vs-legacy
equivalence assertions could be replaced by hand-computed
expressions strictly stronger than legacy-agreement. The B3 sweep
rewrote every such site in one commit cycle (C-B3.7 through
C-B3.12) and deleted the ``apply`` methods.

**β1 interim** (sub-option within Issue #186 B3, considered but
not landed as a final state). Keep
:class:`LinearOperator` inheritance on
:class:`BoundaryTraceLaw` and drop only the abstract
:meth:`apply`. Rank-N composition would then build an
:class:`OperatorSum` of :class:`ScaledOperator`-wrapped *unrealized*
laws, and :func:`realize_recursively` would walk the resulting
Wave-0 composer tree. This β1 form was **algebraically equivalent**
to β2 but conflated two type families: the
:class:`OperatorSum` instance representing a not-yet-realized
descriptor tree was structurally identical to the
:class:`OperatorSum` representing an actual operator composition,
and only runtime inspection of the leaves could tell which was
which. β2 was preferred because the type system can then enforce
"this is a law, that is an operator" by static inspection — a
type checker rejects ``law_tree.apply(...)`` at the linter level
without ever running the code. β2 is more verbose (two new
dataclasses) but is the architecturally-checkable form.

**The vacuum two-paths-divergence is eliminated by design.** Under
B3 + β2 there is no "direct path" any more — every consumer must
realize the law before applying it, and realize-then-apply goes
through the inflow-only-mask path. The Wave-6 snapshot harness's
``vacuum_lebedev17`` case still pins inflow rows only (the
intentional §16A.5 divergence is still there in the algebra), but
no caller can accidentally invoke the pre-§16A.5 zeros-all path
because that path no longer exists in the code.

Empirical justification
-----------------------

The 18-test :mod:`tests.geometry.test_law_composition` suite pins
the descriptor-tree contract (foundation + L1 tests):

* Algebra closure on every dunder for every node-type pairing
  (laws × LawSum × LawScaled).
* :class:`LawScaled` constant folding
  (``2 * (3 * spec) == LawScaled(6, spec)``).
* The walker's exhaustive dispatch on the three node types and
  its :class:`TypeError` on unknown nodes.
* Walker value-correctness against hand-composed expectation
  (``realize_recursively(law_tree).apply(psi)``
  ``== 0.3 * realize(spec).apply(psi) + 0.7 * realize(white).apply(psi)``).
* The absence of ``apply`` on any descriptor-tree node
  (``not hasattr(tree, "apply")``).

The Wave-6 snapshot harness (now realizer-only — the legacy
halves are gone) verifies that the realized output is
bit-identical to the pre-B3 realizer-path output on every
non-mixed case, and identical up to the documented ULP tolerance
on the Marshak mixed case (the operator tree is structurally the
same; only the route to it changed).

Call-site contract
------------------

There is **one** way to call a boundary law's ``apply``:

.. code-block:: python

   from orpheus.geometry.boundary import ReflectiveBoundary
   from orpheus.sn.boundary.realizer import (
       SNBoundaryRealizer, SNMethodSpace,
   )

   law = ReflectiveBoundary(axis="x", albedo=0.5)
   ms = SNMethodSpace.minimal(quad)   # or .for_face(...) in production
   op = SNBoundaryRealizer().realize(law, ms)
   psi_in = op.apply(psi_out)

For descriptor-tree composition:

.. code-block:: python

   from orpheus.geometry.boundary import realize_recursively

   tree = 0.3 * ReflectiveBoundary(axis="x") + 0.7 * WhiteBoundary(
       axis="x", outward_sign=+1,
   )
   op_tree = realize_recursively(tree, ms, SNBoundaryRealizer())
   psi_in = op_tree.apply(psi_out)

No call site routes through a putative ``law.apply(psi)`` — that
method does not exist.


.. _bc-numerical-evidence:

Numerical evidence
==================

The Wave 6 snapshot harness
(:mod:`tests.geometry.test_bc_equivalence_snapshot`) is the
load-bearing equivalence pin for the realizer-vs-legacy migration.
Eight snapshot cases pin pre-refactor outputs as ``.npz`` files;
post-refactor outputs are compared at ``nulp ≤ 4`` (or ``nulp = 64``
for cases where the reduction tree intentionally changed). The
cases:

.. list-table:: Wave 6 BC equivalence snapshot cases
   :header-rows: 1
   :widths: 28 28 18 26

   * - Snapshot
     - BC
     - Quadrature
     - Tolerance
   * - ``vacuum_lebedev17``
     - ``VacuumInflow()``
     - Lebedev order 17
     - inflow-rows-only
       (intentional §16A.5 semantic divergence on outflow)
   * - ``albedo_05_lebedev17``
     - ``AlbedoBoundary(0.5)``
     - Lebedev order 17
     - ``nulp ≤ 4``
   * - ``specular_x_lebedev17``
     - ``ReflectiveBoundary(axis="x", albedo=1.0)``
     - Lebedev order 17
     - bit-identical (α=1 fast path)
   * - ``specular_y_partial_07_LS6``
     - ``ReflectiveBoundary(axis="y", albedo=0.7)``
     - LevelSymmetricSN(6)
     - ``nulp ≤ 4``
   * - ``white_xmax_LS4``
     - ``WhiteBoundary(axis="x", outward_sign=+1, albedo=1.0)``
     - LevelSymmetricSN(4)
     - bit-identical (α=1 fast path)
   * - ``white_xmin_partial_03_GL``
     - ``WhiteBoundary(axis="x", outward_sign=-1, albedo=0.3)``
     - GaussLegendre1D(8)
     - ``nulp ≤ 4``
   * - ``periodic_lebedev17``
     - ``PeriodicBoundary()``
     - Lebedev order 17
     - bit-identical (angular identity body)
   * - ``mixed_30spec_70white_LS4``
     - ``0.3 * spec + 0.7 * white`` (Wave-0 algebra)
     - LevelSymmetricSN(4)
     - ``nulp ≤ 64``

The 1-ULP tolerance for α=1 specular / white reflects the
α=1.0 fast path returning the bare primitive (no
``ScaledOperator`` wrap → same FP reduction tree as the legacy
body). The 4-ULP tolerance for partial-albedo cases reflects the
``ScaledOperator(α, primitive).apply(x) = α * primitive.apply(x)``
multiplication being the single new floating-point operation per
ordinate (one extra rounding step). The 64-ULP tolerance for the
Marshak case reflects the
:class:`OperatorSum` reduction order (the legacy ``MixedBoundaryOperator``
loop summed left-to-right via accumulation; the post-Wave-11
:class:`OperatorSum` binary tree adds in a different but
mathematically equivalent order — see the
``vv-principles`` skill's "bit-identity vs principled-equivalence"
discussion for the general framing of this kind of FP-non-
associativity drift).

The snapshot ``.npz`` files live at
``tests/geometry/snapshots/bc_equivalence_*.npz`` and are committed
to the repository. The generator script
``tests/geometry/_generate_bc_equivalence_snapshots.py`` regenerates
them on the legacy commit; the comparison test enforces
equivalence on every subsequent commit. Wave 11 (``MixedBoundaryOperator``
removal) regenerated the mixed snapshot — the new payload was
**bit-identical** to the pre-Wave-11 payload because the realizer's
deleted internal mixed-BC path was already composing via
``OperatorSum``-of-``ScaledOperator``.


.. _bc-two-bc-applies-per-matvec:

Two BC apply calls per curvilinear matvec
=========================================

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
  outer-face flux.  For vacuum BC the
  :class:`IncomingOrdinateMaskTensor` zeroes the inflow ordinate,
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
calls :meth:`bc.apply` twice and the test must locate Call #2 (the
§16A.3 call) unambiguously.

The Phase D test
:func:`tests.sn.sweep.core.test_phase_c_gates.test_bc_trace_contract_capture_and_compare_sphere`
(parametrised over ``vacuum`` and ``reflective``):

#. Monkey-patches ``sn_mesh.bc["xmax"].apply`` (the outer radial
   face — a sphere's ``"outer"`` endpoint renders as ``"xmax"``)
   with a recorder wrapper that appends every input array passed to
   it during one matvec call.
#. Independently reconstructs the WDD-propagated outflow trace via
   a reference implementation
   (``_outflow_at_boundary_for_sphere(sn_mesh, sig_t, psi_input)``).
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
``vacuum`` case is the load-bearing check because Call #1 produces
non-trivial output under vacuum (the
:class:`IncomingOrdinateMaskTensor` zeroes inflow ordinates, so the
extracted ``bc_outer_value`` is zero — but the **input** to Call #1
is the outer cell-centre value which is **not** zero on a non-trivial
:math:`\psi`).  Locating Call #2 unambiguously requires content
matching, not just shape matching.


.. _bc-three-bc-applies-per-sweep-iteration:

BC applies in the SI sweep path
===============================

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
:mod:`tests.sn.sweep.core.test_sweep_vs_apply_consistency`
(57 foundation tests) pins that the sweep-path's
``bc_outer_value`` extraction matches the apply-path's Phase D
Call #1 result on every test configuration — the structural
invariant that the two paths' Carlson seeds agree on matching
inputs subsumes the BC-apply-input pinning.


.. _bc-curvilinear-realizer-unification:

Curvilinear realizer unification
================================

The pre-cleanup architecture carried a **Cartesian / curvilinear
split** at :meth:`SNMesh._resolve_one`: the slab and 2-D Cartesian
paths constructed a trace space (then named ``InflowTraceSpace``,
unified into :class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace`
post-#188 and made geometry-blind in C5.3) and routed the BC through
:class:`SNBoundaryRealizer`, while spherical and cylindrical ``Mesh1D``
bypassed the realizer entirely because that factory
(``from_mesh_and_quadrature``, since C5.3
:meth:`AngularTraceSpace.from_quadrature_and_layout
<orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.from_quadrature_and_layout>`)
raised :class:`NotImplementedError` on those coord systems. The bypass
wrapped the bare law instance in
:class:`_BoundBoundaryOperator(law, quadrature=self.quad)`, a
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
cylindrical :class:`Mesh2D` (which has no SN sweep in ORPHEUS
today) continues to raise :class:`NotImplementedError`.

Issue #188 / C188.3 then collapsed :meth:`SNMesh._resolve_one`
to a single path: every supported mesh (1-D Cartesian / spherical
/ cylindrical + 2-D Cartesian) builds an
:class:`SNMethodSpace.for_face` and routes through
:class:`SNBoundaryRealizer`. Issue #176 / C176.1 then dropped the
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
*how the set of resolved laws is keyed and stored on the* ``SNMesh``.
Pre-C4 that storage was a hand-listed per-geometry construction with
named attributes — ``bc_xmin`` / ``bc_xmax`` / ``bc_ymin`` /
``bc_ymax`` (2-D), ``bc_left`` / ``bc_right`` aliases (1-D), plus a
pair of degenerate ``bc_ymin`` / ``bc_ymax`` placeholders on a slab
that no production code ever read. Three separate hand-lists carried
the same ``(axis, endpoint) → "{axis}{min|max}"`` knowledge, and a
fourth hand-list mapped a face name back to a reflection axis. C4
(part of the N-D layout campaign, Issue #220) collapses all four to
**one crosswalk function and one dict-comprehension loop**, keyed by
the same :class:`~orpheus.transport.mesh.axis.FaceLabel` inventory the trace
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
     - :attr:`SNMesh.boundary_face_layout`
       (:class:`~orpheus.numerics.face_layout.FaceLayout`)
     - The flat-buffer descriptor: which faces exist, each face's
       per-face shape ``(N, ng, *codim-1 cells)``, and the offsets
       that pack them into one backing array.
   * - **Trace space**
     - :attr:`SNMesh._trace`
       (:class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace`)
     - The inner-product geometry on the boundary: per-face
       inflow / outflow ordinate masks over the signed
       :math:`\Omega\cdot\hat n` it carries (``trace.layout.faces``
       reproduces the same names).
   * - **BC dict**
     - :attr:`SNMesh.bc`
       (``dict[str, _BoundBoundaryOperator]``)
     - The resolved boundary operator per face — the realized
       1-arg law that maps an outgoing face trace to its incoming
       partner.

Pre-C4 each of these grew its key set from its own per-geometry
hand-list. The crosswalk knowledge — "axis 0 is ``x``; the
``min`` / ``max`` endpoints suffix the axis name; a solid radial
axis's single ``outer`` endpoint renders as the ``max`` face"
— was duplicated at the layout builder, the BC resolver, and the
reflective-axis dispatch. C4 lifts that knowledge into **one**
single-sourced rendering on the structural face key:

.. code-block:: python

   # orpheus/sn/axis.py
   AXIS_NAMES = ("x", "y", "z")
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

:attr:`FaceLabel.face_name <orpheus.transport.mesh.axis.FaceLabel.face_name>`
is THE rendering of the structural identity ``(axis_index,
endpoint)`` into the ``"{axis}{min|max}"`` string world. Both
producers — :attr:`SNMesh.boundary_face_layout` and
:meth:`SNMesh.realize_boundary_law` — call it, so a key drift between the
face layout and the BC dict is **unrepresentable by construction**:
they cannot disagree because they read the same function over the
same :func:`~orpheus.transport.mesh.axis.face_labels` inventory.

.. note::

   ``AXIS_NAMES`` moved **down** from
   :mod:`orpheus.sn.loss_representation.sweep_graph` to :mod:`orpheus.transport.mesh.axis` in C4 —
   to the bottom of the SN dependency graph, next to the axis
   primitives it names. ``sweep_graph`` re-exported it only outward;
   ``sweep_schedule`` and ``loss_representation`` now import it
   downward. This puts the single source of the axis↔name crosswalk
   in the same module as :class:`~orpheus.transport.mesh.axis.FaceLabel`, the
   walk's in/outflow-face derivation, and the schedule's
   outgoing-face derivation — no consumer hand-lists
   ``("x", ...), ("y", ...)`` pairs any longer.

The ``"outer" → "max"`` convention and fail-loud
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A solid sphere or cylinder is a :class:`~orpheus.transport.mesh.axis.RadialAxisMesh`
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
:class:`ValueError`. An :class:`~orpheus.transport.mesh.axis.AxisMesh` exposes
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
:attr:`SNMesh.axes` is therefore the root of every boundary-keyed
structure:

.. code-block:: text

   SNMesh.axes  (tuple[Axis1D, ...])
        │
        │  face_labels(axes) — one FaceLabel per (axis, endpoint),
        │  iterated axis-ascending then endpoint-in-axis-order
        ▼
   SNMesh.face_labels  (tuple[FaceLabel, ...])
        │
        ├─── boundary_face_layout : one slot per label,
        │       named  label.face_name,  shaped (N, ng, *face_shape(label))
        │
        └─── bc : one entry per label,
                keyed  label.face_name,  resolved from
                axes[label.axis_index].bc[label.endpoint]

The BC **declaration** source for each face is the per-axis
inventory ``axes[label.axis_index].bc[label.endpoint]`` — the
*same* axes that ``face_labels`` derives the labels from. The face
inventory **is** the BC inventory: a face that exists has exactly
one declaration; a face that does not (the pole) has no label and no
entry. The resolution loop is one comprehension:

.. code-block:: python

   # orpheus/sn/geometry.py — _resolve_bcs (post-C4)
   default = BC("reflective")
   self.bc: dict[str, _BoundBoundaryOperator] = {
       label.face_name: self._resolve_one(
           self.axes[label.axis_index].bc[label.endpoint] or default,
           label,
       )
       for label in self.face_labels
   }

``None`` on an axis defaults to ``BC("reflective")`` (the
infinite-lattice / eigenvalue convention). Each declaration is
realized by :meth:`SNMesh._resolve_one`, whose realizer plumbing
(registry → :meth:`SNMethodSpace.for_face` →
:class:`SNBoundaryRealizer` → :class:`_BoundBoundaryOperator`) is
**unchanged** from the pre-C4 path — C4 changed only the *keying and
storage*, not the *realization*. Hence the resolved operators are
bit-identical objects to the pre-C4 ones (see
:ref:`bc-face-name-carve-verification`).

The :attr:`boundary_face_layout` producer is the dual comprehension:

.. code-block:: python

   # orpheus/sn/geometry.py — boundary_face_layout (post-C4)
   return FaceLayout.from_named_shapes([
       (label.face_name, (N, self.ng, *self.face_shape(label)))
       for label in self.face_labels
   ])

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

A :class:`~orpheus.transport.mesh.axis.RadialAxisMesh` has
``endpoints = ("outer",)`` — exactly one BC-bearing endpoint. A
solid sphere or cylinder therefore has a ``bc`` dict with **exactly
one entry** (``"xmax"``) and **no pole entry**. The geometric pole
at :math:`r=0` is the angular closure's regularity condition (the
:math:`1-\mu^2` redistribution coefficient vanishes at
:math:`\mu=\pm 1`; the inward sweep seeds from a moment-folded
source at :math:`\mu=-1` — see
:ref:`sn-pole-angular-closure-protocol` in
:doc:`/theory/methods/sn/curvilinear_one_group`), not a BC trace law.

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
     - ``sn_mesh.bc_xmin`` → ``None``
     - ``sn_mesh.bc["xmin"]`` → :class:`KeyError`
   * - Failure mode of a buggy consumer
     - silent ``None.apply`` → ``AttributeError`` deep
       in a sweep, or a guard that *should* exist but
       doesn't
     - immediate ``KeyError`` at the access site — the
       illegal access is unrepresentable rather than null

This is the illegal-states-unrepresentable principle (Pattern 4)
applied to the dict: a pole-BC is not "a BC that is ``None``", it is
"not a face at all". The :func:`~orpheus.transport.mesh.axis.face_labels`
inventory simply does not emit a label for the pole, so neither
producer writes a slot or entry for it.

.. _bc-face-name-latent-d3-bug:

The latent d=3 axis-dispatch bug, closed by construction
--------------------------------------------------------

Before C4, :meth:`SNMesh._resolve_one` derived a reflective law's
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

C4 derives the axis from the label's own
``AXIS_NAMES[label.axis_index]``, so the reflection partner is
correct at **any** dimension by construction::

    law = ReflectiveBoundary(axis=AXIS_NAMES[label.axis_index], albedo=1.0)

This is the boundary-resolution sibling of the C3.6 finding that a
z-face never sheds in the in-plane projection — both are latent
3-D correctness traps closed *before* ``Mesh3D`` exists, so that the
N-D layout campaign reaches C5 with the boundary keying already
3-D-correct.

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
   is documented here, not in ``error_catalog.md``.

Why string-keyed, not FaceLabel-keyed
-------------------------------------

Issue #220 allowed either a ``dict[str, ...]`` keyed by face name or
a ``dict[FaceLabel, ...]`` keyed by the structural label directly.
C4 chose **string-keyed**, isomorphic to
:attr:`FaceLayout.faces <orpheus.numerics.face_layout.FaceLayout>`,
for one reason: **every consumer iterates ``trace.layout.faces``
(strings)**. The within-group operator's
:meth:`SNBoundaryOperator._face_laws` and the schedule's
``_reflective_faces`` both walk the trace layout's face-name strings
and index the BC by that string. A ``FaceLabel``-keyed dict would
force a reverse ``name → label`` lookup at *every* consumer, re-deriving
the very crosswalk C4 single-sources.

:class:`~orpheus.transport.mesh.axis.FaceLabel` remains the **structural source**
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
   pinned by the foundation gate ``tests/sn/primitives/test_face_name_crosswalk.py``
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
  aliases (1-D). Consumers now key into :attr:`SNMesh.bc` by face
  name. Accessing a retired attribute is an :class:`AttributeError`.
* **The degenerate 1-D y-face placeholders.** Pre-C4, a slab
  :class:`SNMesh` carried a pair of realized no-op
  ``ReflectiveBoundary(axis="y")`` operators at ``bc_ymin`` /
  ``bc_ymax``, routed through :meth:`SNMethodSpace.minimal` so
  cross-dimensional code could read them without coord-system
  gating. **No production code ever read them**: a 1-D mesh's
  ``trace.layout.faces`` is ``("xmin", "xmax")``, so the generic
  consumers (which iterate the trace layout) never asked for a
  y-face. The placeholders were a uniformity affordance with no
  consumer — exactly the kind of dead realized state the
  face-labels-derived dict makes unrepresentable: a slab has no
  y-axis in its :attr:`~SNMesh.axes` tuple, so
  :func:`~orpheus.transport.mesh.axis.face_labels` emits no y-label, so
  :attr:`bc` has no y-entry, so ``slab.bc["ymin"]`` is a
  :class:`KeyError`. (Pre-C4 design rationale for *why the
  placeholders were once safe* is preserved in the
  :ref:`bc-curvilinear-realizer-unification` history above; C4
  removes the need for the rationale by removing the placeholders.)
* **The hand-listed reflective-axis dispatch**
  (``"y" if face in (...) else "x"``) — see
  :ref:`bc-face-name-latent-d3-bug`.

Consumers migrated in C4: :meth:`SNBoundaryOperator._face_laws`
(within-group operator) and ``sweep_schedule._reflective_faces``
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
  construction; :meth:`_resolve_one`'s realizer plumbing (registry
  → :meth:`SNMethodSpace.for_face` → :class:`SNBoundaryRealizer` →
  :class:`_BoundBoundaryOperator`) is unchanged. The resolved
  operators are the same objects as before, so every solver test
  that exercises them inherits its prior verification. The affine
  ``sha256`` goldens stayed byte-identical; the broad sweep /
  operators / primitives / solve suite is green.
* **L0 crosswalk pins**
  (:mod:`tests.sn.primitives.test_face_name_crosswalk`,
  foundation-tagged). An exhaustive **hand-transcribed**
  ``(axis, endpoint) → face-name`` table for :math:`d \in \{1,2,3\}`
  (mirror-not-import, so the test is not tautological against the
  production derivation it verifies), the d=3 ``z``-face admission
  (the crosswalk is a pure function, so the 3-axis rendering is
  verifiable **now** with no ``Mesh3D``), and both fail-loud
  negatives (a non-canonical endpoint → :class:`ValueError`; an
  axis beyond the named inventory → :class:`IndexError`).
* **L0 bc-dict / face-layout inventory pins**
  (:mod:`tests.sn.operators.test_snmesh_realizer_wiring`,
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
  inheritance from :class:`BoundaryTraceLaw`, and formalised the
  descriptor-tree composition algebra via the new
  :class:`LawSum` / :class:`LawScaled` types. The architectural
  sequence is therefore Issue #188 → #176 → #186: each step
  unblocked the next.

The :class:`_BoundBoundaryOperator` shim survives because the
``kind``-string tag is load-bearing for the BC-resolution
diagnostic and several ``sn_mesh.bc["xmin"] ==
"vacuum"``-style test sites (the dict-keyed spelling since C4 / #220;
this was ``sn_mesh.bc_left == "vacuum"`` pre-C4); the dual-mode
bound-quadrature
backing is gone (#176), and the ``*_extra, **_kw`` swallow is
gone (#186 / C-B3.4). Every supported mesh consumes a strict
1-arg :class:`LinearOperator` produced by
:class:`SNBoundaryRealizer` for single BCs, or by
:func:`realize_recursively` for rank-N descriptor trees.

Plan documents:

* ``.claude/plans/transient-giggling-cake.md`` — the foundational
  12-wave BC trace-law refactor plan (Waves 0–12 close-out
  documented at :ref:`theory-boundary-conditions`).
* ``.claude/plans/curvilinear-realizer-and-2arg-cleanup.md`` —
  the #188 + #176 cleanup plan (Option-A landing).
* ``.claude/plans/bc-trace-law-descriptor-cleanup.md`` — the
  Issue #186 B3 + β2 cleanup plan (descriptor-model landing).

Grand Report v3 §16A.3 (the three-layer architecture) is now
**enforced by the type system** — descriptors have no ``apply``,
operators do. Grand Report v3 §16A.5 (the trace-correct vacuum
representation) is uniform across coord systems and the legacy
zeros-all path no longer exists.


.. _sn-axis-primary-c5:

The axis-primary inversion and 3-D admission
============================================

C4 (:ref:`bc-face-name-carve`) made the *boundary keying*
dimension-agnostic. C5 makes the **whole mesh** dimension-agnostic and
then admits the first 3-axis Cartesian :class:`SNMesh` — *without* a
``Mesh3D`` dataclass. The design fork (resolved by the user,
2026-06-11) is **axis-native**: a 3-D problem enters ORPHEUS only
through :meth:`SNMesh.from_axes` with a 3-tuple of
:class:`~orpheus.transport.mesh.axis.AxisMesh`. :class:`~orpheus.geometry.mesh.Mesh1D`
and :class:`~orpheus.geometry.mesh.Mesh2D` stay the :math:`d \le 2`
user-facing surface, bit-identical to before
(``sha256`` affine goldens unchanged, no regeneration). A ``Mesh3D``
would have had **exactly one consumer** (SN — ``cp`` / ``mc`` / ``moc`` /
``diffusion`` consume zero ``Mesh2D``); the "Unify after two instances"
discipline forbids minting a base type for a single consumer.

The campaign's keystone insight, surfaced by the C5 elegance audit, was
that the d=3 admission could not be a clean *extension* until a
pre-existing **data-flow inversion** in the constructor was repaired. C5
is therefore sequenced *clean before extend*: C5.1–C5.4 invert and
de-phantom the mesh layer, and only then C5.5 admits d=3 as a
one-line gate removal.

.. _sn-c5-lossy-roundtrip:

Pre-C5: the lossy axes → mesh → axes round-trip
-----------------------------------------------

The SN phase space factors as a tensor product of per-axis 1-D meshes
(grand report §15.1); the natural primary representation of an
:class:`SNMesh` is therefore its **axes tuple**
:attr:`SNMesh.axes`. Pre-C5.1 the constructor did not treat it that
way. :meth:`SNMesh.from_axes` *synthesized a legacy*
:class:`~orpheus.geometry.mesh.Mesh1D` / :class:`~orpheus.geometry.mesh.Mesh2D`
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
   :class:`~orpheus.transport.mesh.axis.AxisMesh` carries user-overridable
   ``label_low`` / ``label_high`` fields (a slab user may name them
   ``"left"`` / ``"right"``). The legacy mesh has no slot for those
   labels, so the round-trip dropped them and the re-derived axes came
   back with default labels — a silent desync of exactly the kind C4's
   :attr:`FaceLabel.face_name <orpheus.transport.mesh.axis.FaceLabel.face_name>`
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

   SNMesh(mesh, ...)   ──►  axes = axes_from_legacy_mesh(mesh)   (convert ONCE,
                                                                  at the inbound
                                                                  boundary)
   from_axes(axes, ...) ─►  axes  (stored verbatim — the caller's tuple)
                                            │
                                            ▼
                                       _init_core(axes, ...)

The legacy ``SNMesh(mesh, ...)`` surface converts via
``axes_from_legacy_mesh`` **once**, at the inbound boundary
(parse-don't-validate); :meth:`SNMesh.from_axes` stores the caller's
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
   ``tests/sn/primitives/test_axis_native_construction.py``
   (``test_d2_metadata_byte_identical_axis_vs_legacy`` /
   ``test_1d_slab_metadata_byte_identical_axis_vs_legacy``). A carve
   representational identity, not a solver claim.
.. vv-status: sn-axis-widths documented

This is **bitwise identical** to the legacy per-dataclass spellings it
replaces — :attr:`Mesh1D.widths <orpheus.geometry.mesh.Mesh1D>`,
:attr:`Mesh2D.dx <orpheus.geometry.mesh.Mesh2D>`, and
:attr:`Mesh2D.dy <orpheus.geometry.mesh.Mesh2D>` are each
``np.diff(edges)`` over the same edge arrays (``mesh.py:287`` /
``:567`` / ``:572``), so the carve produces the same floating-point
bytes. The whole-mesh coordinate system is likewise derived from the
per-axis coordinates by a new pure primitive
:func:`~orpheus.transport.mesh.axis.coord_system` (a multi-axis mesh must be
all-Cartesian); the constructor's reduced-operator dispatch and the
pole-closure default now read the **axis-derived** :attr:`SNMesh.coord`,
not ``mesh.coord``.

After C5.1, :attr:`SNMesh.mesh` is **inbound provenance only** — it
records *which legacy mesh the caller passed, if any*. It is ``None``
when the mesh was built from axes at :math:`d \ge 3` (no legacy mesh
exists to record). A handful of :math:`d \le 2` consumers (the 1-D
reduced streaming constructors, the trace build, realizer metadata)
still read ``self.mesh`` at C5.1; each dissolves across C5.2–C5.5 as
its datum is repointed to an axis-native source. ``legacy_mesh_from_axes``
narrows from a round-trip *source* to a :math:`d \le 2` **adapter**
synthesis for those remaining consumers.

Custom endpoint labels now fail loud (C4 doctrine)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

With the axes stored verbatim, a custom endpoint label survives
construction — and therefore reaches the
:attr:`FaceLabel.face_name <orpheus.transport.mesh.axis.FaceLabel.face_name>`
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

With per-axis widths and the rank-generic :attr:`SNMesh.spatial_shape`
now native, the legacy phantom-bearing metadata retires. Every spatial
read becomes rank-honest:

.. list-table:: The retired spellings and their rank-honest replacements
   :header-rows: 1
   :widths: 24 30 46

   * - Retired spelling
     - The phantom it carried
     - Replacement
   * - ``SNMesh.ny`` / ``SNMesh.dy``
     - At :math:`d = 1` these **lied** — ``ny`` returned a phantom
       ``1`` and ``dy`` a phantom ``[1.0]`` (the Issue #214 phantom
       class), and at :math:`d \ge 3` they underspecify the mesh.
     - :attr:`SNMesh.spatial_shape` (the per-axis cell counts) and
       :attr:`SNMesh.axis_widths` (per-axis widths). ``AttributeError``
       on the retired names.
   * - ``SNMesh.dx``
     - A duplicate spelling of the per-axis widths.
     - :attr:`SNMesh.axis_widths` — promoted from the private
       ``_axis_widths`` to **the** single public spelling of per-axis
       cell widths.
   * - ``SNMesh.nx``
     - (kept) Documented :attr:`spatial_shape[0] <SNMesh.spatial_shape>`
       sugar — honest at any :math:`d`, with a broad legitimate 1-D
       consumer base.
     - unchanged.

The phantom ``ny`` / ``dy`` at :math:`d = 1` is the **same defect
class** the N-D layout campaign closed at C2 / Issue #214: a trailing
singleton that masquerades as a real axis. The C5.2 retirement removes
the masquerade at the metadata source.

The two production ``dr`` consumers (the
:mod:`~orpheus.sn.loss_representation` 1-D bare sweep and the
:mod:`~orpheus.sn.sweep.pole_angular_closure` Carlson preamble) repoint
from ``.dx`` to :attr:`SNMesh.axis_widths`. The
field / cross-section / scattering read-through chains collapse to the
rank-generic :attr:`spatial_shape <SNMesh.spatial_shape>`:

* :class:`~orpheus.transport.fields.angular_flux.AngularFlux` (and the
  ``BulkField`` base) **retire** their ``nx`` / ``ny`` read-throughs.
  This is a live :math:`d = 3` landmine, not cosmetic: an
  ``(nx, ny)``-keyed field read **silently truncates** a 3-D tensor to
  its first two axes (a ``vv-principles`` Mode-2 / Mode-5 class
  index error that the degenerate :math:`d \le 2` test never reaches).
* :class:`~orpheus.transport.mesh.material_xs_field.MaterialXSField` and
  :class:`~orpheus.transport.operators.scattering.ScatteringOperator` collapse their
  ``nx`` / ``ny`` reads to **one** rank-generic
  :attr:`spatial_shape <SNMesh.spatial_shape>` read-through each.

Finally, a new :attr:`SNMesh.volume_measure` property gives the
SN-side ``keff`` rate consumers (the production / absorption rates in
:mod:`~orpheus.sn.solver`) a native source: they read it instead of
reaching through ``sn_mesh.mesh.volume_measure``. While the
:math:`d \le 2` adapter is present it delegates to the dataclass
measure (bit-identical, including the ``precomputed_volumes`` hatch);
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
derives the table from :data:`~orpheus.numerics.face_layout.AXIS_NAMES`
so every axis-aligned face (all **six** at :math:`d = 3`) is present by
construction: face ``"{axis}min"`` has outward normal
:math:`-\hat e_{\text{axis}}`, face ``"{axis}max"`` has :math:`+\hat
e_{\text{axis}}`. The :math:`\Omega\cdot\hat n` row for ``zmax`` is
then exactly :math:`+\mu_z` and for ``zmin`` exactly :math:`-\mu_z`.

To share the ``"{axis}{min|max}"`` crosswalk without an ``sn``-ward
import from the ``numerics`` layer, ``AXIS_NAMES`` **moved down** to
:mod:`orpheus.numerics.face_layout` — the home of
:class:`~orpheus.numerics.face_layout.FaceLayout`, keeper of the face
string-name world. :mod:`orpheus.transport.mesh.axis` re-exports it, so SN consumers
are unchanged; the trace space (a ``numerics`` leaf) now reads it
without depending on ``sn``.

**The trace factory is geometry-blind.** The mesh parameter on the old
``from_mesh_and_quadrature`` factory was **gate-only** — its single use
was an ``isinstance`` check that refused a curvilinear ``Mesh2D``. That
refusal is **unreachable**: a curvilinear ``Mesh2D`` cannot become an
:class:`SNMesh` in the first place (2-D cylindrical SN has no sweep), so
no such mesh ever reached the factory. The ``isinstance`` check carried
no data the factory used. C5.3 therefore renames the factory to
:meth:`AngularTraceSpace.from_quadrature_and_layout
<orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.from_quadrature_and_layout>`
and **retires the dead mesh parameter** (aggressive retirement; callers
and the bare-constructor error message migrated). Every datum now comes
from the quadrature (the :math:`\mu_x` / :math:`\mu_y` / :math:`\mu_z`
cosines) and the layout's face names (the axis-aligned normals implied
by the ``"{axis}{min|max}"`` convention).

With the gate gone, :meth:`SNMesh.realize_boundary_law` builds the trace
**unconditionally** (the pre-C5.3 ``isinstance`` gate excluded only the
unconstructible 2-D cylindrical mesh), and :attr:`SNMesh.angular_trace` is typed
and documented as **always non-None**.
:meth:`SNMethodSpace.for_face <orpheus.sn.mesh.method_space.SNMethodSpace.for_face>`'s
``mesh`` parameter becomes **optional metadata** — nothing in the
realizer chain reads it, and an axis-native :class:`SNMesh` passes
``None``.

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
``sn_mesh.reduced is None``:

* the **moment-windowing** gate (:meth:`_maybe_window
  <orpheus.sn.solver>`), which decides whether the SI iterate is held
  as compact harmonic moments rather than the full angular flux; and
* the **Gauss–Seidel resolvent** selector (``_select_si_resolvent``),
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
   * - Moment windowing (``_maybe_window``)
     - ``reduced is None``
     - ``is_cartesian and ndim == 2`` — the genuine
       windowing-eligibility condition (the 2-D moment kernel's exact
       domain).
   * - Boundary-G-S (``_select_si_resolvent``)
     - ``reduced is None``
     - ``is_cartesian and not is_1d`` — multi-D Cartesian.

The G-S resolvent's old ``"2-D Cartesian ONLY"`` docstring was **stale
Phase-3 narration**: :attr:`SweepSchedule.gauss_seidel
<orpheus.sn.loss_representation.sweep_schedule>` and the scheduled sweep
(``_sweep_scheduled``) have been **d-generic since C3**, so the
resolvent is constructible at :math:`d = 3`. The narration is corrected
in C5.4; the actual :math:`d = 3` boundary-G-S *fixed-point invariance*
(that G-S and Jacobi converge to the **same** flux) is **value-gated**
by the C5.5 Mode-9 mixed-BC box (:ref:`sn-c5-value-gates`) before any
:math:`d = 3` G-S solve is trusted — the Mode-9 discipline made
operative: never gate a splitting's FP-invariance on a degenerate box.

The one ``reduced is not None`` branch that **stays** is the 1-D
sweep-cache: there the predicate keys on the *availability of the data
it reads* (the reduced-operator cache), not on dimensionality. That is
not a proxy — it is the genuine guard.

.. _sn-c5-3d-admission:

3-D Cartesian admission: the axes tuple is the only 3-D entry
-------------------------------------------------------------

After the C5.1–C5.4 cleanup, the :math:`d = 3` admission is an
**extension, not a new arm**. A 3-axis Cartesian :class:`SNMesh` now
constructs and **solves** through the same generic body as
:math:`d \le 2`, **mesh-adapter-free from birth** (``self.mesh is
None``) on the d-generic
:class:`~orpheus.sn.loss_representation.FullFieldWavefront` spine.

* **The gate retires.** :meth:`SNMesh.from_axes` drops the
  ``d \ge 3`` admission guard; :math:`d \le 2` still synthesizes the
  legacy adapter for its remaining consumers.
* **Axis-native arms.** The cell-volume array is the iterated outer
  product of the per-axis widths,
  :math:`V[i,j,k] = \mathrm{d}x_i\,\mathrm{d}y_j\,\mathrm{d}z_k`; the
  :attr:`volume_measure <SNMesh.volume_measure>` is the rank-:math:`d`
  meshgrid-of-centers
  :class:`~orpheus.numerics.measure.DiscreteMeasure` (the natural
  rank-:math:`d` generalization of the ``Mesh2D`` analogue).
* **Entry surface.** :func:`~orpheus.sn.solver.solve_sn` and
  :func:`~orpheus.sn.solver.solve_sn_fixed_source` accept the **axes
  tuple** — the *only* 3-D entry — through one inbound seam
  (``_as_sn_mesh``). A new ``mat_map`` keyword is the axes-entry
  material channel (it raises if combined with a legacy mesh, which
  carries its own material map). Default-BC semantics are handled per
  surface (``_apply_default_bcs`` accepts both declaration styles —
  per-face dataclass fields *or* per-endpoint axis slots — with the
  same all-or-nothing semantics).

.. note::

   **Two default-BC conventions, by design.** The *solver* entry
   defaults un-declared faces to **vacuum** (the fixed-source
   convention — an un-specified boundary leaks); a freshly constructed
   :class:`SNMesh` with no BC declarations defaults to **reflective**
   (the infinite-lattice / eigenvalue convention — see
   :ref:`bc-face-name-carve`). The d=3 admission preserves **both**
   conventions on their respective surfaces; the value gates below
   exercise the reflective (eigenvalue) convention for the headline
   :math:`k_\infty` identity and a mixed convention for the Mode-9 box.

.. _sn-c5-value-gates:

Numerical evidence — the d=3 value gates
----------------------------------------

C5.5's admission is gated by four value tests
(:mod:`tests.sn.solve.test_d3_admission`), **all driven through the
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
   * - **Mode-9 G-S ≡ Jacobi FP-invariance**
     - L2 (integration)
     - Boundary-Gauss–Seidel and Jacobi converge to the **same** d=3
       fixed point on a box that **breaks every degenerate
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

Separately, the **multi-D adjoint** (``loss_action_transpose`` raises
:class:`NotImplementedError` at any multi-D) is a **pre-existing
deferral** orthogonal to C5 (G-adjoint campaign territory), noted here
only so the :math:`d = 3` capability map is complete.

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
re-attempt them:

1. **Single ``BoundaryOperator`` ABC carrying both law and
   realizer responsibilities.** This is the pre-refactor shape.
   Rejected because the law / realizer split is the architectural
   point of the refactor: laws are method-agnostic, realizers are
   method-specific. Keeping them in one class would force every
   law to know about the SN sweep's inflow-mask plumbing, which is
   precisely what Wave 0 / Wave 1 / Wave 2 was supposed to abstract
   away.
2. **Dedicated ``MixedBoundaryOperator`` class.** Pre-Wave-11
   shape. Rejected (deleted Wave 11) — the Wave-0 algebra dunders
   on every :class:`BoundaryTraceLaw` already produce
   :class:`OperatorSum` shapes; the dedicated class added a
   second realization path with no semantic difference. See
   :ref:`bc-rank-n-algebra`.
3. **Shared ``BoundaryRealizerBase`` ABC for cross-method
   realizers.** Considered Wave 5, when only one functional
   realizer existed. Rejected per the "Unify after two instances"
   architectural discipline: building the abstraction on a single
   instance would force a particular shape on every future method
   based on SN's dispatch idiom (``isinstance``). Vindicated at
   method #2: the diffusion realizer (#290 P3) chose a DIFFERENT
   shape (law → albedo scalar → structure-keyed collapse, not an
   isinstance ladder over per-law primitives), and the structural
   :class:`~orpheus.geometry.boundary.BoundaryRealizer` Protocol
   remains the only shared contract — no ABC was ever needed.
4. **Adding ``face`` to ``VacuumInflow``'s constructor for
   semantic correctness on the standalone-apply path.** Option (b)
   in :ref:`bc-vacuum-semantic-correction`. Rejected because it
   would have inflated every test signature for one wave to fix a
   path that the refactor was retiring anyway. The transitional
   legacy-vacuum body is the right cost.
5. **Auto-importing every cross-method realizer at
   :mod:`orpheus.geometry.boundary` import time.** Considered in
   the registry era to make ``BoundaryRealizerRegistry.get("MoC")``
   work without the caller having to ``import orpheus.moc`` first.
   Rejected because :mod:`orpheus.sn` is a heavy module that every
   consumer of the boundary package would then pay for. #290 P7b
   made the whole question moot by dissolving the registry: each
   method-mesh imports its own realizer explicitly, so there is no
   name-lookup to keep populated and no import-side-effect timing
   to defend.
6. **Cartesian-vs-curvilinear bypass in
   :meth:`SNMesh._resolve_one` + dual-mode shim.** Pre Issue #188
   shape: curvilinear ``Mesh1D`` bypassed the realizer and wrapped
   the bare 2-arg law in
   ``_BoundBoundaryOperator(law, quadrature=self.quad)``, while
   Cartesian routed through the realizer with
   ``_BoundBoundaryOperator(realized)``. **Retired Issue #188 +
   #176**; documented here because the seductive trap is to
   "preserve flexibility" by keeping the dual mode after the
   curvilinear deferral lifts. The right move is to delete the
   bypass and consolidate on one path — see
   :ref:`bc-curvilinear-realizer-unification`.
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
8. **Calling ``apply`` on a raw BC descriptor.** Under the
   pre-#186 contract this either worked (with surprising
   semantics — see Option A entry above) or raised
   :class:`BoundaryError`. Under post-#186 it's a **static type
   error** — :class:`BoundaryTraceLaw` has no :meth:`apply`
   method on the class, and neither do :class:`LawSum` /
   :class:`LawScaled`. The correct contract is
   ``SNBoundaryRealizer().realize(law, ms).apply(psi)`` for a
   single BC, or ``realize_recursively(tree, ms).apply(psi)`` for
   a descriptor tree. The realizer is the **sole** bridge; the
   §16A.3 three-layer split is enforced by the type system.
9. **In-tree Wave-0 operator algebra over unrealized
   :class:`BoundaryTraceLaw` instances (β1 form).** Considered as
   the rank-N composition mechanism during Issue #186 B3
   exploration. Rejected in favour of the separate-type-family
   approach (β2 / :class:`LawSum` / :class:`LawScaled`). β1
   produced :class:`OperatorSum` trees whose leaves were laws,
   not operators — the type checker could not distinguish a
   not-yet-realized "operator" from a real operator, and calling
   :meth:`apply` on the β1 tree raised at the leaf realization
   step. β2 makes the law-vs-operator distinction inspectable
   statically: :class:`LawSum` has no :meth:`apply` method, full
   stop. See :ref:`bc-rank-n-algebra` for the detailed
   comparison.


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
  :math:`R = \sum_\alpha G_\alpha \otimes A_\alpha` that this page
  refines into the affine form :eq:`affine-bc-form`.
* :ref:`operator-algebra` for the Wave-0 primitives the realized
  BCs decompose into.
* The V&V error catalog in the ``vv-principles`` skill
  (``.claude/skills/vv-principles/error_catalog.md``) carries
  the ERR-040..ERR-047 entries in canonical form.
