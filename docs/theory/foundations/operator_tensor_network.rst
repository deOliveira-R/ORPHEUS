.. _operator-tensor-network:

=====================================================
Tensor-Network Decomposition of S\ :sub:`N` Operators
=====================================================

.. contents:: Contents
   :local:
   :depth: 2


.. Machine header — the ``nexus-meta`` schema for this page (PROVISIONAL).
.. This page was extracted verbatim from ``operator_algebra.rst``; the schema
.. is provisional pending a full re-audit of the split corpus (#231).

.. dropdown:: Machine header — ``nexus-meta`` schema (PROVISIONAL)
   :color: muted

   .. code-block:: yaml

      module: transport
      concept: operator_tensor_network
      role: "tensor-network / factored-shape decomposition of the S_N transport operators"
      depends_on: [operator_algebra]
      status: "extracted from operator_algebra.rst; content verbatim, provisional header"


This page develops the **factored / tensor-network shape decomposition**
of the S\ :sub:`N` transport operators: *which* algebraic shape each
operator leaf actually takes once it is lifted out of a procedural,
single-axis numpy body into the operator-algebra types
(:class:`~orpheus.numerics.operator.TensorProductOperator`,
:class:`~orpheus.numerics.operator.SumOfTensorProductsOperator`,
:class:`~orpheus.numerics.operator.OperatorSum`) developed in
:doc:`/theory/foundations/operator_algebra`. The headline result is that
**no single uniform tensor-product factorisation across space, angle and
energy** exists: the shipped state is *five algebraically distinct shapes*,
chosen per operator by what the underlying physics actually couples — and
the streaming operator, whose in-sweep recurrence couples the axes
sequentially, resists a clean tensor product altogether.

- :ref:`wave-t-shape-table` — the per-operator shape catalogue and the
  MA-Q1 master condition that decides when a tensor product is admissible.
- :ref:`wave-t-streaming-deep-dive` — streaming's in-sweep WDD recurrence
  and why it is not tensor-separable.
- :ref:`wave-t-orchestrated-apply` — the single bidirectional pass, its
  design rationale, and the retired per-direction split (#238).
- :ref:`wave-t-curvilinear-deep-dive` — the curvilinear angular
  redistribution (Morel–Montry) thread.


.. _wave-t-tensor-network:

Tensor-Network Decomposition of SN Operators
============================================

Wave T (May 2026, commits ``fa13e78`` / ``0b2848b`` / ``9f85c5d`` /
``03bcdba`` / ``cb18fdb`` / ``c55b505`` / ``90e7d4e``) lifted the four
SN operator leaves — boundary realizers, fission, scattering,
streaming — from procedural single-axis numpy bodies into the
operator-algebra types documented in
:doc:`/theory/foundations/operator_algebra`
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


Key Facts
---------

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
     - :class:`OperatorProduct`
       :math:`R \circ \Lambda_{\ell\ge 1} \circ M`
     - the moment-space integral kernel
       ``frame.conjugate(Λ)``
       (:attr:`ScatteringOperator.kernel`), projecting to harmonic
       moments, applying per-ℓ transfer, reconstructing per-ordinate
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
       :doc:`/theory/methods/diffusion_1d`), inverted by an explicit
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

Why the split was retired — the orphan-smell one level down
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
   :doc:`/theory/methods/diffusion_1d` — per the architecture decided on
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
:doc:`/theory/methods/sn/curvilinear_one_group`, Step 2–3,
Eq. :eq:`balance-general`), it
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
:class:`orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase`
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
experiments (recorded during the pre-reorganisation ``orpheus/sn/operator.py``
work, since split into ``orpheus/sn/operators/``). That rewire requires the BC
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
     ``tests/sn/verification/analytical/test_kinf_homogeneous.py``. This is the
     eigenvalue reference — MMS does NOT prove eigenvalues per
     the ``vv-principles`` skill §"What each pillar
     proves".

   - **MMS pillar**: P1 anisotropic manufactured-source convergence
     at ``tests/sn/verification/mms/test_mms_aniso.py``,
     ``tests/sn/verification/mms/test_curvilinear_aniso_convergence.py``,
     ``tests/sn/verification/mms/test_mms_heterogeneous.py``, and
     ``tests/sn/verification/mms/test_mms_2d.py``. The MMS source is structurally
     independent of the operator-algebra path (derived by SymPy in
     ``orpheus/derivations/continuous/mms/sn.py``); it catches flux-shape and
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
:class:`OperatorSum` over per-ℓ ``_PerLegendreOrderScattering``
bespoke leaves. The §15.2 *form* is preserved at the summation level
(one summand per Legendre order); the per-summand decomposition into
:math:`R_\ell \circ \Lambda_\ell \circ M_\ell` is a procedural
composition, not a tensor product. (This per-ℓ ladder was itself
later retired at commit ``93807aa7`` in favour of the single
Funk–Hecke :math:`R \circ \Lambda_{\ell\ge 1} \circ M` moment-space
kernel now carried by :attr:`ScatteringOperator.kernel`.)

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
    :attr:`~orpheus.transport.operators.scattering.ScatteringOperator.kernel`
    property.
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
  - :class:`orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase`
    — the M-M closure data and per-cell algebra primitive.
  - :func:`orpheus.transport.spatial.cell_balance.cell_balance_for_streaming`
    — the three-term cell-balance primitive.
