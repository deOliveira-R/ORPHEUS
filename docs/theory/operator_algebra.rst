.. _operator-algebra:

=============================
Operator Algebra Architecture
=============================

ORPHEUS' transport, eigenvalue, and Krylov solvers all act on a flux
distribution by composing a small set of linear operators — the
streaming + collision operator :math:`L`, the scattering operator
:math:`S`, the fission operator :math:`F`. The
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
  :class:`~orpheus.transport.source_sinks.boundary_source_sink.BoundarySourceSink`
  (Wave O step B.5.2, Issue #208, ``6ef5063``, 2026-06-03), completing
  the boundary half of the operator-output "dimensional-sin" carve (the
  bulk half — ``.apply.bulk`` →
  :class:`~orpheus.transport.source_sinks.angular_source_sink.AngularSourceSink`
  — landed in ``f400743``). The governing principle: an operator output
  is :math:`A\psi` — a **source/sink**, NOT a residual; the residual
  arises ONLY from an explicit :meth:`~orpheus.transport.residuals.boundary_residual.BoundaryResidual.from_balance`
  of the output against a source. The completed boundary role grid
  mirrors the bulk:
  :math:`\texttt{.apply} \to \texttt{BoundarySourceSink}`,
  :math:`\texttt{.solve} \to \texttt{BoundaryFlux}`,
  :math:`\texttt{from\_balance} \to \texttt{BoundaryResidual}`. See
  :ref:`bc-extraction-operator-output-typing`.

- **The interior cell-face angular fluxes are a typed cochain**
  :class:`~orpheus.transport.fields.wavefront_flux.WavefrontFlux`
  (Wave O step #205 Phase 5, Issue #208, ``478723d`` / ``992b0c0`` /
  ``0e3e16c``, 2026-06-04): the 2-D wavefront sweep and matvec no
  longer carry raw ephemeral ``psi_x`` / ``psi_y`` numpy arrays. The
  interior 1-cochain :math:`C^1_{\rm int}` and the boundary trace
  :class:`~orpheus.transport.fields.boundary_flux.BoundaryFlux`
  (the boundary 1-cochain :math:`C^1_\partial`) **biproduct-decompose
  the full face cochain** :math:`C^1 = C^1_{\rm int} \oplus
  C^1_\partial` — the :math:`V_{\rm bulk} \oplus V_{\rm boundary}`
  shape of :eq:`bc-extraction-direct-sum-state` one locus down, at the
  *face* level. The seed/absorb the sweep applies by hand are now the
  typed trace operators :math:`\iota_*` (:meth:`~orpheus.transport.fields.wavefront_flux.WavefrontFlux.seed`)
  and :math:`\iota^*` (:meth:`~orpheus.transport.fields.wavefront_flux.WavefrontFlux.absorb`),
  with the "absorption = identity" fact now the provable biproduct law
  :math:`\iota^* \circ \iota_* = \mathrm{id}`. See
  :ref:`wavefront-flux-cochain`.

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
  :class:`~orpheus.sn.boundary_operator.SNBoundaryOperator` — the same
  operator the 2-D Krylov path uses. The legacy "B1'' face block"
  (never a code symbol; a 1-D boundary-closure *name*) is retired. See
  :ref:`bc-extraction-2d-si-krylov-twin`.

- The transport eigenvalue and fixed-source problems share a single
  operator algebra (Trefethen & Bau 1997, §3.2):

  .. (vv-status rationale) Phase 0 stub label for the transport
     fixed-source operator-algebra form. The verifiable claim — that
     a Phase-1 SN/MoC/CP solver assembles its operator triple
     (L, S, F) and lets ``op = L - S - F`` agree with the legacy
     fixed-source path bit-for-bit — is queued for Issues 9-13.
  .. vv-status: operator-fixed-source documented

  .. math::
     :label: operator-fixed-source

     (L - S - F)\,\psi \;=\; q
     \qquad \text{(fixed source)}

  .. (vv-status rationale) Phase 0 stub label for the transport
     eigenvalue operator-algebra form. Verified by Issue 14's unified
     iteration driver consuming (L - S, F) triples, and by the
     existing eigenvalue test suites once they are migrated to the
     new contract.
  .. vv-status: operator-eigenvalue documented

  .. math::
     :label: operator-eigenvalue

     (L - S)\,\psi \;=\; \tfrac{1}{k}\,F\,\psi
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
  admits monolithic-matrix resolvents that have no :math:`(L,S,F)`
  triple);
  :class:`~orpheus.numerics.iteration.KEigenvalue` is the
  operator-triple realization that *delegates* its loop to it (one
  loop in the codebase, Cardinal Rule 2). See
  :ref:`eigenvalue-posing`.

- The :class:`~orpheus.numerics.operator.LinearOperator` Protocol
  carries one mandatory method (``apply``) and a
  ``capabilities: frozenset[str]`` property advertising which optional
  methods (``solve``, ``apply_transpose``) are functional.

- **SN array-storage convention** for every operator leaf
  (:class:`~orpheus.sn.operator.StreamingOperator`,
  :class:`~orpheus.sn.operator.CollisionOperator`,
  :class:`~orpheus.sn.scattering.ScatteringOperator`,
  :class:`~orpheus.sn.fission.FissionOperator`): the
  ``apply(psi) -> psi'`` contract consumes and returns
  ``psi.shape == (N, ng, nx, ny)`` for angular flux,
  ``phi.shape == (ng, nx, ny)`` for scalar flux.  The canonical
  statement with derivation and migration history lives at
  :ref:`theory-sn-index-convention`.

- The capability set is the **single source of truth** for what an
  operator can do. Composition primitives compute their own capability
  set from constituents; mismatches raise
  :class:`~orpheus.numerics.operator.MissingCapability` at composition
  time, NEVER at call time.

- :func:`~orpheus.numerics.operator.as_scipy_linop` adapts any
  :class:`LinearOperator` to scipy's Krylov solvers (BiCGSTAB, GMRES)
  via :class:`scipy.sparse.linalg.LinearOperator` — preserving the
  capability set as the gate that decides whether ``rmatvec`` is
  exposed.


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
— only that ``solve`` advertises a meaningful approximation. Operators
without an efficient ``solve`` simply do not advertise that
capability, and downstream code declines them at composition time.


Capability set semantics
========================

The capability set is :class:`frozenset[str]` rather than a class
hierarchy because many transport operators have **no** efficient
``solve``:

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
discover the failure mid-iteration. The capability-set contract is
"I can do *these* things" rather than "I can do everything or
fail" — composers veto incompatible compositions at construction
time.

The set is open in the sense that user operators MAY advertise
method-specific capability tags (e.g. a custom preconditioner could
add ``"jacobi-preconditioned"``); composers ignore tags they do not
understand, propagating only the three primitive ones.


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

The closure rules above are enforced at composition time via
:class:`~orpheus.numerics.operator.MissingCapability`; a downstream
consumer that requests ``solve`` on :math:`L = A + 0` cannot silently
hit a meaningless answer because the sum has already vetoed
``solve`` propagation.


Cross-solver consumption (forward reference)
============================================

The Phase 0 algebra is the foundation for the cross-solver migration
sequence: SN, MoC, CP, and diffusion will each expose their natural
operator triple :math:`(L, S, F)` as
:class:`~orpheus.numerics.operator.LinearOperator` instances, and a
unified iteration driver in :mod:`orpheus.numerics.iteration` will
consume those triples without knowing which transport discretisation
they came from. That driver is **to be installed** under the Wave H
refactor (Issue 14). The Phase 0 module ships the algebra; subsequent
phases mount the producers and the consumer onto it.

Today the existing
:func:`build_transport_linear_operator <orpheus.sn.operator.build_transport_linear_operator>`
functions in :mod:`orpheus.sn.operator` continue to wrap their matvec
in :class:`scipy.sparse.linalg.LinearOperator` directly, untouched by
this module. The :func:`as_scipy_linop` adapter exists so that when
those call sites migrate to the ORPHEUS protocol (Issue 15), the
scipy interop continues to work without rewriting BiCGSTAB / GMRES
inner loops.


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
     - This work, Wave 1 (used inside
       :class:`HarmonicMomentProjection`).
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
operator advertises ``CAP_SOLVE`` and its :meth:`solve`
divides by :math:`w_n` along ``axis``. If any weight is zero,
``solve`` is dropped from the capability set at construction time
— a zero weight has no inverse, and the harmful-stub anti-pattern
the capability-set design exists to prevent (a downstream Krylov
consumer silently dividing by zero) is caught upfront.

The implementation does NOT eagerly materialise an
:math:`N \times N` diagonal matrix. The action is a single
broadcast-multiply
(``w.reshape((1, ..., -1, ..., 1)) * x``) so memory cost is
:math:`O(N)` regardless of the input tensor's shape. For the SN
angular axis this matters: an :math:`(N, n_x, n_y, n_g)` field with
:math:`N = O(10^2)` and :math:`n_x \cdot n_y \cdot n_g = O(10^7)`
does NOT need a :math:`(10^7, 10^7)` materialised diagonal.


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
   product reference on small concrete factors, capability-
   intersection, and the algebraic laws below.
.. vv-status: tensor-product-action documented

.. math::
   :label: tensor-product-action

   (A_1 \otimes A_2 \otimes \cdots \otimes A_k)\,x
   \;=\; A_k\bigl(\cdots A_2(A_1\,x) \cdots\bigr).

Because the constituents act on disjoint axes, the order does not
matter — the operators commute on the joint tensor. The
``capabilities`` set is the **intersection** of the
constituents' capabilities: the tensor product can apply iff every
factor can apply, can apply_transpose iff every factor can, can
solve iff every factor can.


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
     - Typed operators; ``axis`` attribute on each factor; capability
       intersection; algebraic laws checked at composition.
     - ``HarmonicMomentProjection(weights=w, Y=Y, L=L) &
       IdentityOperator()`` — the type signal carries the axis
       structure; mismatched axes raise at composition; the
       :meth:`apply` routes through ``np.einsum`` internally with
       the correct subscripts.

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
* The capability set is computed automatically; a missing
  capability raises at composition, not at the first
  :func:`numpy.einsum` call mid-iteration.

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

For per-octant operators (e.g. the SN sweep's :math:`L_{oct}^{-1}`
acting on
:math:`(\text{octant\_ordinates} \times \text{cells} \times
\text{groups})`), the **tensor product** factors the per-axis
structure within an octant while the **direct sum** assembles the
octants into the global angular cubature:

.. math::
   :label: octant-direct-sum-tensor-product

   L^{-1} \;=\; \bigoplus_{\text{oct}}\,
                L_{oct}^{-1}, \qquad
   L_{oct}^{-1} \;=\; \text{(per-axis tensor product within an octant)}.

.. vv-status: octant-direct-sum-tensor-product documented

This is the operator-algebra type signature behind Wave 2 of the
SN performance plan: the 2-D Cartesian sweep iterates
:math:`\bigoplus_{oct}` (4 octants — structural) and per-octant
anti-diagonals (sweep-DAG topology — structural), with no Python
loop over the ordinate axis (which is **internal** to every
:meth:`apply` call within an octant, vectorised by the tensor-
product structure).


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
   :class:`LinearOperatorMixin` inheritance was dropped). The
   table below lists the **realised** output produced by
   :class:`~orpheus.sn.boundary_realizer.SNBoundaryRealizer` — that
   output IS a Wave-0 :class:`LinearOperator`. The law-to-operator
   transition is the §16A.3 three-layer architecture's realiser
   bridge, enforced statically by the type system rather than by
   convention. See :ref:`bc-trace-law-descriptor-model` for the
   design rationale.

Three new primitives in :mod:`orpheus.numerics.operator`
(Wave 0 of the BC refactor) plus one SN-specific primitive in
:mod:`orpheus.sn.angular_operator` (Wave 1) form the realisation
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
   * - :class:`~orpheus.sn.angular_operator.AngularAverageOperator`
     - Cosine-weighted Lambertian average over an outgoing
       hemisphere: scalar-broadcasts to all inflow ordinates.
       SN-specific (depends on
       :class:`~orpheus.sn.quadrature.AngularQuadrature`).
     - :class:`~orpheus.geometry.boundary.WhiteBoundary`
   * - :class:`~orpheus.sn.angular_operator.IncomingSourceOperator`
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
:func:`~orpheus.sn.boundary_realize.realize_recursively` walks the
tree, realises each leaf via
:class:`~orpheus.sn.boundary_realizer.SNBoundaryRealizer`, and
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
  :class:`LinearOperatorMixin` on :class:`BoundaryTraceLaw`, so
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

- **Per-direction streaming split**:
  :attr:`StreamingOperator.M_spatial` is an :class:`OperatorSum` of
  forward (μ_x > 0) and backward (μ_x < 0) sweep summands — NOT a
  monolithic sweep. The split is architecturally load-bearing for
  Wave O typing, adjoint propagation, DSA-class preconditioners,
  cross-method (billiard / trajectory-resolvent) pollination, and
  per-direction debugging slicing.

- **Orchestrated apply (Design B)**:
  :class:`_MSpatialOperatorSum` overrides the default
  :class:`OperatorSum`'s additive apply because the two per-direction
  summands share local state — the forward-sweep outer-face WDD
  outflow seeds the backward sweep's :attr:`bc_outer` boundary
  condition. The naïve sum would cost 1.5× the unified matvec because
  each standalone leaf would re-run the forward sweep to compute its
  own seed. Perf measured at median 1.04× pre-T.4 baseline (commit
  ``90e7d4e``).

- **Hybrid 2-D Cartesian**: T.4 lifted 1-D only. The 2-D Cartesian
  path (:meth:`StreamingOperator._apply_2d_cartesian`) stays
  procedural FD with cell-centre-proxy semantics. A defensive
  source-hash pin (A2D-1) guards against silent author drift on that
  path until a future wave delivers the trace-correct face_view
  refactor.

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
       (single rank-1 TP)
     - ``RankOneOperator(χ, νΣ_f, axis=0) &
       IdentityOperator()``
       (:attr:`FissionOperator.kernel`)
     - Per Grand Report v3 §15 line 2008
       :math:`F = \chi \otimes \nu\Sigma_f`. The
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
       part
       (:attr:`StreamingOperator.M_spatial`)
     - :class:`OperatorSum` of two
       per-direction bespoke leaves
       (subclass
       :class:`_MSpatialOperatorSum`)
     - ``_SpatialSweepDirection(+1) +
       _SpatialSweepDirection(-1)`` orchestrated
       via :meth:`_MSpatialOperatorSum.apply`
     - **MA-Q1 fallback**: the WDD recurrence
       :math:`\psi_{\text{face,out}} = 2\bar\psi
       - \psi_{\text{face,in}}` sequentially
       couples cells along x. The per-direction
       summands expose structure at the type
       level but are NOT clean
       :math:`(D_x \otimes \Omega_x \otimes I_g)`
       3-factor TPs — the sweep operator is the
       leaf factor.
   * - Streaming angular
       redistribution
       (:attr:`StreamingOperator.M_angular_redist`)
     - Bespoke
       :class:`LinearOperator`
       leaf (sphere / cylinder)
       OR :class:`ZeroOperator`
       (slab / 2-D Cartesian)
     - :class:`AngularRedistributionOperator`
       wrapping
       :meth:`PoleAngularClosure.cell_contribution`
     - **MA-Q1 fallback**: the M-M half-grid
       recurrence (Hébert 2009 §3.9.4
       Eqs. 3.432-3.435) sequentially couples
       angular ordinates ``α_{m+1/2}`` from
       ``α_{m-1/2}`` with σ_t-dependent
       absorption coefficients. Not a diagonal
       angular factor; a 3-factor TP wrap would
       false-assert separability the recurrence
       doesn't support.


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

Streaming :math:`M_{\rm spatial}` deep dive — per-direction split
-----------------------------------------------------------------

:attr:`StreamingOperator.M_spatial` is an :class:`OperatorSum` of
**two per-direction-sign summands**:
:class:`_SpatialSweepDirection(+1)` (forward sweep, :math:`\mu_x > 0`)
and :class:`_SpatialSweepDirection(-1)` (backward sweep,
:math:`\mu_x < 0`). The split serves five distinct future consumers,
none of which justifies the split on its own but which together carry
significant architectural payload.

**Why split per direction:**

1. **Wave O typing** — the natural type signal for the per-direction
   summand exposes the BC dependency footprint cleanly: the forward
   summand reads :attr:`bc_inner` (or the symmetry axis for sphere /
   cylinder) and writes the outer-face WDD residual; the backward
   summand reads the outer-face inflow (the BC's
   :meth:`apply(outgoing)` result) and writes the inner-face WDD
   residual. Without the split, every consumer would have to
   re-derive this dependency structure.

2. **Adjoint propagation**. The reverse-direction sweep is the
   structural transpose of the forward sweep with swapped boundary
   conditions: :math:`(M_{x,+})^* = M_{x,-}^{\text{BC-swapped}}`. The
   per-direction split exposes the dual operator's identity at the
   type level, so the adjoint sensitivity path
   (deferred to Phase H) can match per-direction summands without
   mining the internal sweep body.

3. **DSA-class preconditioners**. Traditional diffusion-synthetic
   acceleration schemes split per-direction at the P_1 closure level.
   The per-direction summand is the structural unit the DSA
   preconditioner consumes.

4. **Cross-method pollination**. The billiard solver and the
   trajectory-resolvent reference both expose per-direction sweep
   primitives in their own architectures (a single forward-trajectory
   evaluation in billiard; a single ray-direction Bickley-Naylor pass
   in trajectory-resolvent). Exposing the per-direction summand at
   the SN type level provides a future structural bridge.

5. **Per-direction debugging slicing**. The ERR-006 family
   (curvilinear sweep divergence) and Signature-1 recurrence-coupled
   bugs (catalogued in the ``vv-principles`` skill) manifest
   *per direction* before they manifest in the summed
   matvec. A per-direction property gives test code a clean
   inspection point without re-implementing the sweep body in the
   test.

**WDD coupling along x**. The sequential coupling that breaks SOTP
applies *within* each per-direction summand. The forward sweep is

.. math::
   :label: wdd-forward-recurrence

   \bar\psi_i \;=\;
     \frac{\Delta x_i\,\bar Q_i + |\mu|\,(\psi_{\text{face,in}})_i}
          {\Delta x_i\,\Sigma_t(i) + |\mu|},
   \qquad
   (\psi_{\text{face,out}})_i \;=\;
     2\,\bar\psi_i \;-\; (\psi_{\text{face,in}})_i

with :math:`(\psi_{\text{face,in}})_{i+1} =
(\psi_{\text{face,out}})_i`. The cell-balance algebra at
:func:`orpheus.sn.spatial.cell_balance.cell_balance_for_streaming`
hides this recurrence inside the named denom-numer primitives, but
the recurrence is the load-bearing structure that prevents a clean
:math:`(D_x \otimes \Omega_x \otimes I_g)` 3-factor TP.

.. vv-status: wdd-forward-recurrence documented

**Forward-backward coupling at the outer face**. The two
per-direction summands are independent *in their per-direction
contribution*, but they share local state: the forward sweep's
outer-face WDD outflow IS the input to the backward sweep's outer-face
BC application. In the legacy unified-matvec body, this shared state
was a local variable; under the per-direction operator-algebra split,
it becomes the **named shared state** of the
:class:`_MSpatialOperatorSum` orchestrator.

.. note::

   The forward-outflow-seeds-backward-inflow coupling described in
   this paragraph was the *pre*-O.4a.2 mechanism, where the backward
   sweep's inflow was ``bc_outer.apply(forward_outflow)`` inside one
   matvec. Wave O step O.4a.2 (Issue #208) **deleted that intra-call
   reflective re-apply** — the boundary law :math:`B` is now a sibling
   :math:`-B` operator and the backward sweep reads the *given* outer
   inflow trace directly. The orchestrator's named shared state today
   is the forward outflow as it feeds the **outflow self-consistency
   defect** on the outflow trace row, not a reflected inflow seed.
   See :ref:`bc-extraction`.


.. _wave-t-orchestrated-apply:

Orchestrated :meth:`apply` (Design B)
-------------------------------------

:class:`_MSpatialOperatorSum` is a **subclass** of
:class:`OperatorSum` (so its ``.a`` and ``.b`` attributes carry the
two :class:`_SpatialSweepDirection` summands for type introspection
by Wave O / adjoint / DSA). The subclass **overrides**
:meth:`OperatorSum.apply` because the default

.. math::

   \texttt{OperatorSum.apply}(x) \;=\;
     \texttt{self.a.apply}(x) \;+\; \texttt{self.b.apply}(x)

would cost 1.5× the unified matvec walltime: each standalone
:meth:`_SpatialSweepDirection.apply` invokes
:meth:`_MSpatialOperatorSum._compute_LpC <orpheus.sn.operator._MSpatialOperatorSum._compute_LpC>` (which already runs the
**bidirectional** sweep internally) and masks the opposite-direction
ordinates to zero. Calling the unified matvec twice — once per
direction summand — duplicates the forward sweep cost.

The orchestrator runs the bidirectional sweep **once** via
:meth:`_MSpatialOperatorSum._compute_LpC <orpheus.sn.operator._MSpatialOperatorSum._compute_LpC>` and returns the full
:math:`(L+C)\,\psi` for slab (since slab has no curvilinear
redistribution; see :ref:`wave-t-curvilinear-deep-dive` below for the
curvilinear subtraction). The standalone
:meth:`_SpatialSweepDirection.apply` is preserved as a slow fallback
for testing, Wave-O adjoint inspection, and per-direction debugging
slicing.

**Why this matters architecturally**. The forward-sweep outer-face
WDD outflow that was a hidden local variable in the legacy unified
matvec body is now the **named shared state** of the orchestrator
(:meth:`_MSpatialOperatorSum._compute_LpC <orpheus.sn.operator._MSpatialOperatorSum._compute_LpC>`).
Pattern 6 (single source of
truth) of the project's ``coding-elegance`` skill requires that hidden
coupling points become named. Wave T does not refactor the sweep
body; it lifts the hidden coupling into a named property at the
operator-algebra level.

The full bidirectional matvec is mathematically equivalent to
:math:`M_{x,+}\,\psi + M_{x,-}\,\psi`, and the
:meth:`_MSpatialOperatorSum.apply` orchestrator returns the same
value bit-exact (Route A — preserving the unified matvec's reduction
order). The standalone per-direction summands return masked outputs
that, when summed, equal the unified matvec output at
FP-non-associativity ULP.


.. _wave-t-curvilinear-deep-dive:

Curvilinear :math:`M_{\rm angular\_redist}` — bespoke leaf
----------------------------------------------------------

For sphere / cylinder geometries, :attr:`StreamingOperator.M_angular_redist`
returns an :class:`AngularRedistributionOperator` — a bespoke
:class:`LinearOperator` leaf wrapping the M-M (Morel-Montry) half-grid
recurrence. The leaf consumes the per-cell M-M algebra at
:meth:`PoleAngularClosure.cell_contribution` (Pattern 6 — single
source of truth for the M-M coefficients).

**Why a leaf and not a tensor product**. Per Hébert 2009 §3.9.4,
Eqs. 3.432-3.435, the M-M closure produces an angular recurrence

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
respects the disjoint-axes contract.

.. vv-status: mm-half-grid-recurrence documented

A 3-factor TP wrap of the form ``leaf & I_x & I_g`` would
**false-assert separability** the recurrence doesn't support:
:meth:`TensorProductOperator.assert_separable` would erroneously
pass for an operator whose leaf factor secretly couples
:math:`(N, n_x)` jointly through the sequential per-level recurrence.
Per the ``coding-elegance`` skill's Pattern 4 (make illegal states
unrepresentable), the converse holds: do not represent states
(separability) that aren't actually legal.

**Per-cell algebra**. The leaf's :meth:`apply` walks every
:math:`(p,\,i)` pair (μ-level × spatial cell) and calls
:meth:`PoleAngularClosure.cell_contribution`. The cell-balance
algebra at
:func:`orpheus.sn.spatial.cell_balance.cell_balance_for_streaming`
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

The :math:`M_{\rm spatial}` summand carries the streaming +
collision-share contribution; the :math:`M_{\rm angular\_redist}`
leaf carries

.. math::

   m_{\rm angular\_redist} \;=\;
     \frac{1}{V_i}\bigl[
        {\rm angular\_denom\_term} \cdot \psi_{\rm cell}
      - {\rm angular\_numer\_upstream}
     \bigr]

with :math:`{\rm angular\_denom\_term} = (\Delta A / w)\,c_{\rm out}`
and :math:`{\rm angular\_numer\_upstream} = (\Delta A / w)\,c_{\rm in}\,
\psi_{m-1/2,\,i,\,g}` per the M-M closure (see
:class:`orpheus.sn.spatial.pole_angular_closure.PoleAngularClosure`
for the closure data).

.. vv-status: wave-t-cell-balance-three-terms documented

**Boundary residual semantics**. :attr:`M_angular_redist` writes
**only to bulk**; angular redistribution is an interior-cell operation
that doesn't traverse the spatial boundary. Output ``boundary`` is the
zero :class:`~orpheus.transport.fields.boundary_flux.BoundaryFlux`.
Only :class:`_SpatialSweepDirection` writes non-trivial face
residuals.


Curvilinear :math:`M_{\rm spatial}` via subtraction
---------------------------------------------------

For curvilinear (sphere / cylinder), :meth:`_MSpatialOperatorSum.apply`
computes :math:`M_{\rm spatial}` by **subtracting**
:math:`M_{\rm angular\_redist}` from the unified
:math:`(L+C)\,\psi`:

.. math::
   :label: wave-t-mspat-curvilinear-subtraction

   M_{\rm spatial}\,\psi \;=\; (L+C)\,\psi \;-\;
                               M_{\rm angular\_redist}\,\psi
   \qquad (\text{curvilinear})

.. vv-status: wave-t-mspat-curvilinear-subtraction documented

This subtraction introduces a minor architectural smell: the curvilinear
:math:`M_{\rm spatial}` depends on :math:`M_{\rm angular\_redist}` for
its definition. The alternative — re-implementing the spatial-only
sweep without M-M coupling at the leaf level — would duplicate the
per-level :func:`dag_walk` + Carlson coupled-pole structure that
already lives in :meth:`_MSpatialOperatorSum._compute_decomposition <orpheus.sn.operator._MSpatialOperatorSum._compute_decomposition>`. The
subtraction is the cleanest path that preserves Pattern 6 (single
source of truth).

The smell is bounded by an **algebra-decomposition invariant test**:

.. math::

   (L+C)\,\psi \;\equiv\;
     M_{\rm spatial}\,\psi \;+\; M_{\rm angular\_redist}\,\psi

at principled-equivalence ULP per
the ``vv-principles`` skill (~16×ULP for the
1-D curvilinear matvec). If the subtraction ever falsely cancels (or
adds) a term, this test fires.


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
regression test that records the source-code signature of
:meth:`_apply_2d_cartesian` and asserts it remains unchanged. The pin
catches a future author who accidentally modernises the 2-D path
while editing nearby code, before the change ships and breaks the
known-good cell-centre-proxy semantics.


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
     scattering kernel, T.4 :attr:`M_spatial`) — see the
     "out of scope" note in
     ``.claude/plans/wave_t_tensor_network.md`` §6 T.5.

   - :math:`(L+C)\,\psi \equiv M_{\rm spatial}\,\psi +
     M_{\rm angular\_redist}\,\psi` at Route A bit-identity (slab) or
     ~16×ULP principled equivalence (curvilinear).

   - :math:`(L+C).{\rm solve}(q)` bit-identical pre/post-Wave-T,
     verifying the WDD sweep procedural inverse was NOT touched
     (the :class:`InvertibleOperator.solve` body is the procedural
     algorithm at :func:`orpheus.sn.sweep.transport_sweep`).

5. **Performance regression gate**. The 1-D slab Krylov benchmark
   measured median 1.04× pre-T.4 baseline (under the 5% threshold).
   The :func:`cached_property` decorators on :attr:`M_spatial` and
   :attr:`M_angular_redist` ensure construction happens once per
   :class:`StreamingOperator` lifetime, not per :meth:`apply`.


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
  - :class:`orpheus.sn.boundary_realizer.SNBoundaryRealizer` —
    the BC realizer dispatching the T.1 lifts.
  - :class:`orpheus.sn.fission.FissionOperator` and its
    :attr:`~orpheus.sn.fission.FissionOperator.kernel` property.
  - :class:`orpheus.sn.scattering.ScatteringOperator` and its
    :attr:`~orpheus.sn.scattering.ScatteringOperator.kernel` /
    :attr:`~orpheus.sn.scattering.ScatteringOperator.kernel_summands`
    properties; the bespoke
    :class:`orpheus.sn.scattering._PerLegendreOrderScattering` leaf.
  - :class:`orpheus.sn.operator.StreamingOperator` and its
    :attr:`~orpheus.sn.operator.StreamingOperator.M_spatial` /
    :attr:`~orpheus.sn.operator.StreamingOperator.M_angular_redist`
    properties; the bespoke
    :class:`orpheus.sn.operator._SpatialSweepDirection`,
    :class:`orpheus.sn.operator._MSpatialOperatorSum`, and
    :class:`orpheus.sn.operator.AngularRedistributionOperator`
    leaves.
  - :meth:`orpheus.sn.operator._MSpatialOperatorSum._compute_LpC`
    (single-emission hot path) and
    :meth:`orpheus.sn.operator._MSpatialOperatorSum._compute_decomposition`
    (dual-emission ``(M_spatial, M_angular_redist)`` split) — the
    unified 1-D matvec body. The module-level
    ``_transport_operator_matvec_unified`` helper was **deleted** at
    Wave T step T.5.2 (commit ``ad813fd``); its body was inlined into
    these two private orchestrator methods. The public surface is the
    :meth:`StreamingOperator.apply` /
    :attr:`StreamingOperator.M_spatial` boundary.
  - :class:`orpheus.sn.spatial.pole_angular_closure.PoleAngularClosure`
    — the M-M closure data and per-cell algebra primitive.
  - :func:`orpheus.sn.spatial.cell_balance.cell_balance_for_streaming`
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
(:class:`~orpheus.transport.fields.boundary_flux.BoundaryFlux`),
partitioned per face by the sign of :math:`\Omega\cdot\hat n` into the
inflow (:math:`\Omega\cdot\hat n < 0`) and outflow
(:math:`\Omega\cdot\hat n > 0`) ordinate slots (the
:class:`~orpheus.numerics.spaces.trace_space.TraceSpace` selectors
:meth:`~orpheus.numerics.spaces.trace_space.TraceSpace.inflow_indices_for_face`
/ :meth:`~orpheus.numerics.spaces.trace_space.TraceSpace.outflow_indices_for_face`,
single source of truth — see :ref:`trace-spaces-doc`).

.. vv-status: bc-extraction-direct-sum-state documented

This is the realisation, for the boundary block, of the Wave T
prediction (:ref:`wave-t-tensor-network`): "Wave O typing must accept
non-SOTP summands." :class:`~orpheus.sn.boundary_operator.SNBoundaryOperator`
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
       (:class:`~orpheus.sn.operator.StreamingOperator`,
       via :class:`~orpheus.sn.operator.InvertibleOperator` ``L+C``)
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
       (:class:`~orpheus.sn.operator.CollisionOperator`,
       :class:`~orpheus.sn.scattering.ScatteringOperator`,
       :class:`~orpheus.sn.fission.FissionOperator`)
     - ``BULK``
     - Bulk → bulk only.
     - :math:`A_{bb}` only; the boundary block is zero.
   * - :math:`B`
       (:class:`~orpheus.sn.boundary_operator.SNBoundaryOperator`)
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
an :class:`~orpheus.sn.angular_operator.AngularAverageOperator` for
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
:meth:`SNBoundaryOperator._apply_faces <orpheus.sn.boundary_operator.SNBoundaryOperator>`
**projects** the emission onto the codomain row — ``apply`` writes the
``inflow_indices_for_face`` slots; ``apply_transpose`` writes the
``outflow_indices_for_face`` slots. *Empirically confirmed before the
fix*: the outflow slots carried nonzero :math:`R\cdot\psi.\text{inflow}`.

**3. The bare sweep seeds inflow from** :math:`\text{rhs.boundary}`,
**not** :math:`\text{initial\_guess.boundary}`.
Under the extraction the WDD sweep ``(L+C).solve`` is **bare** (see
:ref:`bare-sweep-extraction` in :doc:`discrete_ordinates`): it reads
the seeded inflow trace directly instead of re-applying ``bc``.
:meth:`InvertibleOperator._solve_timed_full_field <orpheus.sn.operator.InvertibleOperator._solve_timed_full_field>`
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
:math:`(L, S, F)` operator *triple*. Wave O step O.2a generalised both
:class:`~orpheus.numerics.iteration.SourceIteration` and
:class:`~orpheus.numerics.iteration.KrylovAcceleration` to the
**variadic** shape :math:`\text{Driver}(L_{\rm resolvent},\,*\text{gains})`:
one invertible resolvent :math:`L` plus a homogeneous bag of lagged
coupling operators :math:`g_i`. The two consume the gains identically —

.. math::
   :label: bc-extraction-variadic-matvec

   \text{matvec} \;=\; L.\text{apply} \;-\; \sum_i g_i.\text{apply}
   \,,\qquad
   \text{rhs} \;=\; q_{\rm ext} \;+\; \sum_i g_i.\text{apply}\,.

.. vv-status: bc-extraction-variadic-matvec documented

The driver is now **problem-type-agnostic**: it sees only the
invertible resolvent it must invert and a bag of operators it must lag.
*Which* leaves are gains is a **posing-layer** decision, not an
iteration-layer one (see :ref:`eigenvalue-posing`) — the gains are
exactly the posing's coupling terms.

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

assembled once by the single-source-of-truth helper
:func:`_within_group_triple <orpheus.sn.solver._within_group_triple>`,
which returns ``(L+C, S, B)`` — the invertible resolvent
(:class:`~orpheus.sn.operator.InvertibleOperator`, ``.solve`` = the WDD
sweep), the bulk scattering gain
(:class:`~orpheus.sn.scattering.ScatteringOperator`,
:attr:`block_role <orpheus.numerics.operator.BlockRole>` ``BULK``), and
the boundary reflection gain
(:class:`~orpheus.sn.boundary_operator.SNBoundaryOperator`,
``block_role`` ``BOUNDARY``). The :math:`B\,\psi.\text{outflow}` term
lands on :math:`\text{rhs.boundary}`, which the bare ``(L+C).solve``
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
:eq:`bc-extraction-variadic-matvec` shows. The fixed :math:`(L, S, F)`
triple gave :math:`S` and :math:`F` named slots the *resolvent layer*
never uses — it was encoding a posing-layer role assignment (which
operator is loss, which is the eigen-operator) at the iteration layer,
where it does not belong. Collapsing the triple to a homogeneous
:math:`*\text{gains}` bag moves the role distinction back to the
posing layer (its proper home) and lets a fourth gain (a future
:math:`B`-trace term, an :math:`\alpha`-time term) slot in as a data
addition rather than a new named slot. Existing positional
:math:`(L, S, F)` callers stay source-compatible — ``gains = (S, F)``.

**Why** :math:`B` **is a SEPARATE gain, not folded into** :math:`S`.
Two structural reasons forbid the old fold:

#. **The adjoint metric lives on the trace.** :math:`B` lives on the
   boundary trace (:attr:`domain <orpheus.sn.boundary_operator.SNBoundaryOperator.domain>`
   ``= sn_mesh.trace``), and the cosine-weighted
   :math:`|\Omega\cdot\hat n|\,w` adjoint metric (Wave O step O.2 — the
   codomain inner product of :math:`L`'s boundary-trace block) lives on
   that **trace** domain, not the bulk. Folding :math:`B` into the
   bulk :math:`S` would erase the trace typing the future adjoint
   ``.H`` needs.
#. :math:`B` **cannot join the** :math:`L+C` **preconditioner.**
   :class:`~orpheus.numerics.operator.OperatorSum` does **not** carry
   ``CAP_SOLVE`` (it defines only ``apply`` / ``apply_transpose``), so
   an :math:`L + C - B` sum would strip the ``.solve`` (sweep) that the
   SI step and the Krylov preconditioner depend on. :math:`B` must stay
   a *gain* (lagged, applied) — never a summand of the resolvent.

The old fold type-checked **only because**
:attr:`ScatteringOperator.domain` is ``None`` (it inherits the
:class:`~orpheus.numerics.operator.LinearOperatorMixin` default; bulk
operators type via ``block_role``, not a bulk function space). The
:class:`~orpheus.numerics.operator.OperatorSum` domain-compatibility
check fires only when both operands declare non-``None`` domains that
differ; with :math:`S` untagged the check skipped, so the
trace-typed :math:`B` summed silently with the bulk :math:`S`. Giving
:math:`S` a non-``None`` bulk :class:`~orpheus.numerics.operator.FunctionSpace`
domain — a defensive Pattern-4 typing-completion seam — would make
``OperatorSum`` *reject* a re-introduced :math:`S + B` fold at
construction. That tripwire is **not load-bearing now** (the fold is
gone; nothing sums :math:`S` with :math:`B` again), so the seam is
deferred (see :ref:`bc-extraction-scope-future`).

.. note::

   The drivers' :class:`~orpheus.numerics.iteration.KrylovAcceleration`
   matvec :eq:`bc-extraction-variadic-matvec` and the
   :class:`~orpheus.numerics.iteration.SourceIteration` rhs are now the
   *honest* :math:`(L+C-S-B)\,\psi` and :math:`q_{\rm ext}+S\psi+B\psi`
   — the reassociation :math:`L-(S+B)\to(L-S)-B` is documented as a
   **principled-equivalence** change in
   :ref:`bc-extraction-numerical-evidence` (criterion 3 of the
   ``vv-principles`` bit-identity-vs-principled-equivalence gate), not
   a bug.


.. _bc-extraction-two-routes:

The two :math:`-B` delivery routes
----------------------------------

The same :math:`-B` coupling reaches the sweep two ways, both calling
the **identical**
:class:`~orpheus.sn.boundary_operator.SNBoundaryOperator` (single
source of truth, Cardinal Rule 2):

.. list-table:: The two delivery routes for :math:`-B`
   :header-rows: 1
   :widths: 22 40 38

   * - Route
     - Mechanism
     - Used by
   * - **Variadic gain**
       (:func:`_within_group_triple <orpheus.sn.solver._within_group_triple>`
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
from :meth:`SNBoundaryOperator._reflect_trace <orpheus.sn.boundary_operator.SNBoundaryOperator>`
(:ref:`bc-extraction-reflect-trace`).


.. _bc-extraction-reflect-trace:

The trace-only :math:`A_{ss}` leaf — :meth:`reflect_into_inflow`
----------------------------------------------------------------

:math:`B` is the :math:`A_{ss}` block :math:`V_{\rm outflow} \to
V_{\rm inflow}`: it maps the *outflow* trace to the *inflow* trace.
Both delivery routes ultimately need the same per-face action — apply
each face's realized law (the specular
:class:`~orpheus.numerics.operator.PermutationOperator` for reflective,
:class:`~orpheus.sn.angular_operator.AngularAverageOperator` for
white, zero for vacuum) and project onto the inflow row. To guarantee
they cannot drift, that action is the single
:meth:`SNBoundaryOperator._reflect_trace <orpheus.sn.boundary_operator.SNBoundaryOperator>`
core, and both the full-field forward action
:meth:`B.apply <orpheus.sn.boundary_operator.SNBoundaryOperator.apply>`
and the new trace-only leaf
:meth:`B.reflect_into_inflow <orpheus.sn.boundary_operator.SNBoundaryOperator.reflect_into_inflow>`
route through it (Wave O step O.2a, commit ``8563f4b``).

The leaf exists because the direct helper does not need a full field.
:meth:`B.apply` operates on a :class:`~orpheus.transport.timed_full_field.TimedFullField`
(zero bulk, trace populated) — the bulk is only a carrier to reach the
:math:`A_{ss}` boundary block. The pre-extraction direct helper
fabricated a throwaway zero-bulk field purely to call ``B.apply`` and
then discarded the (zero) bulk output.
:meth:`reflect_into_inflow <orpheus.sn.boundary_operator.SNBoundaryOperator.reflect_into_inflow>`
takes a bare :class:`~orpheus.transport.fields.boundary_flux.BoundaryFlux`
trace and returns the boundary-only
:class:`~orpheus.transport.source_sinks.boundary_source_sink.BoundarySourceSink`
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
:meth:`_reflect_trace <orpheus.sn.boundary_operator.SNBoundaryOperator>`
projects the forward action onto ``inflow_indices_for_face`` and the
Euclidean transpose onto ``outflow_indices_for_face``.


.. _bc-extraction-scope:

Scope — both 1-D and 2-D are now bare (O.4b complete)
-----------------------------------------------------

O.4a.2 made the **1-D** sweep bare (slab / sphere / cylinder). Step
**O.4b** then made the **2-D Cartesian wavefront sweep bare as well**
(both :func:`~orpheus.sn.sweep._sweep_2d_wavefront` and the 2-D matvec
:meth:`StreamingOperator._apply_2d_cartesian <orpheus.sn.operator.StreamingOperator>`):
the intra-sweep ``bc.apply`` is **gone** for every geometry. The
octant-incoming face edge is seeded from the *given* inflow trace and
the reflective coupling :math:`\psi.\text{inflow} = B\,\psi.\text{outflow}`
is delivered by the sibling :math:`-B`
(:class:`~orpheus.sn.boundary_operator.SNBoundaryOperator`) for the
2-D trace exactly as for the 1-D trace. The 2-D matvec emits the
boundary block as an active-trace residual (outflow slots carry the
self-consistency defect ``streamed − ψ.outflow``; inflow slots carry
the identity ``ψ.inflow``), wired into the composed Krylov matvec as
the boundary gain :math:`B` of
:func:`_within_group_triple <orpheus.sn.solver._within_group_triple>`.
The interior face fluxes the bare 2-D sweep + matvec propagate are now
the typed cochain :class:`~orpheus.transport.fields.wavefront_flux.WavefrontFlux`
(see :ref:`wavefront-flux-cochain`).

The dispatch is still guarded by a **single predicate** so the two
geometry paths cannot drift: ``sn_mesh.reduced is not None`` is the
**same** predicate :func:`~orpheus.sn.sweep.transport_sweep` uses to
select the 1-D scan vs the 2-D wavefront body, and the **same**
predicate the direct-helper guards
(:func:`_solve_fixed_source_si <orpheus.sn.solver._solve_fixed_source_si>`,
:func:`solve_sn <orpheus.sn.solver.solve_sn>`) check before calling
:func:`_reflect_outflow_into_inflow <orpheus.sn.solver._reflect_outflow_into_inflow>`.
Both branches are now bare-sweep + sibling :math:`-B`; the predicate
selects the *fold shape* (1-D parallel-prefix scan vs 2-D wavefront
DAG), **not** a bare-vs-bc-in-sweep distinction.


.. _bc-extraction-scope-future:

Deferred typing-completion seam — :attr:`ScatteringOperator.domain`
-------------------------------------------------------------------

One Wave-O typing-completion remains a documented seam, not a feature:
minting a non-``None`` bulk
:class:`~orpheus.numerics.operator.FunctionSpace` domain for
:class:`~orpheus.sn.scattering.ScatteringOperator` (and the other bulk
leaves). Today bulk operators type via ``block_role``, not a domain
space, so :attr:`ScatteringOperator.domain` is ``None``. Giving it a
bulk :math:`V_{\rm bulk}` domain would let
:class:`~orpheus.numerics.operator.OperatorSum` **reject** any attempt
to re-introduce the :math:`S + B` fold — the domain-compatibility check
would throw ``IncompatibleOperatorComposition`` because :math:`S`
(bulk space) and :math:`B` (trace space) live on different function
spaces (a defensive Pattern-4 illegal-states-unrepresentable typing).

This is **not load-bearing now**: with the variadic driver
(:ref:`bc-extraction-variadic-driver`) the fold is gone and nothing
sums :math:`S` with :math:`B` again, so there is nothing for the
tripwire to catch. The seam is recorded so a future typing wave lands
it as a pure addition rather than discovering the need under a
regression.


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
bulk + :meth:`BoundarySourceSink.zeros_on <orpheus.transport.fields._bases.BoundaryField.zeros_on>`
boundary inside a
:class:`~orpheus.transport.timed_full_field.TimedFullField`), the same
loss decomposition (the resolvent :math:`L + C` from
:class:`~orpheus.sn.operator.StreamingOperator` +
:class:`~orpheus.sn.operator.CollisionOperator`, plus the scattering
gain :math:`S` and the boundary reflection gain :math:`B` from
:func:`_within_group_triple <orpheus.sn.solver._within_group_triple>`;
zero within-group fission), and the
same :meth:`integrate_angular <orpheus.transport.fields.angular_flux.AngularFlux.integrate_angular>`
angular reduction. Neither driver carries any geometry dependence.

The reflective coupling reaches both drivers on the **bare** 2-D
wavefront sweep through the sibling :math:`-B` (the **variadic-gain**
route of :ref:`bc-extraction-two-routes`), never through an in-sweep
``bc.apply``. The :class:`~orpheus.sn.boundary_operator.SNBoundaryOperator`
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
   :class:`~orpheus.transport.fields.boundary_flux.BoundaryFlux` +
   :class:`~orpheus.sn.boundary_operator.SNBoundaryOperator`
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
:class:`~orpheus.transport.fields.boundary_flux.BoundaryFlux` to
:class:`~orpheus.transport.source_sinks.boundary_source_sink.BoundarySourceSink`
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
       (minted, no consumer until O.2)
   * - **boundary** (:math:`V_{\rm inflow} \oplus V_{\rm outflow}`)
     - :class:`~orpheus.transport.source_sinks.boundary_source_sink.BoundarySourceSink`
       (``6ef5063``, B.5.2)
     - :class:`~orpheus.transport.fields.boundary_flux.BoundaryFlux`
     - :class:`~orpheus.transport.residuals.boundary_residual.BoundaryResidual`
       (minted, no consumer until O.2)

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
:meth:`BoundaryResidual.from_balance(lhs, rhs) <orpheus.transport.residuals.boundary_residual.BoundaryResidual.from_balance>`
of the affine boundary balance
:math:`r_\Gamma = \gamma_-\psi - (R\,G\,\gamma_+\psi + q)`. The GMRES
*flat* residual :math:`b - A\psi` is formed internally on the **raveled
vector** (via :meth:`TimedFullField.to_flat <orpheus.transport.timed_full_field.TimedFullField.to_flat>`)
and is **never typed as a field** — so :class:`BoundaryResidual` has no
operator-output consumer; its first consumer is the honest
:math:`L+C-S-F-B` driver of Wave O step **O.2** (see
:ref:`bc-extraction-operator-output-o2`). Until then
:class:`BoundaryResidual` and its bulk sibling
:class:`~orpheus.transport.residuals.angular_residual.AngularResidual`
remain **minted, units-tagged, and consumerless** role-grid
completions.


The two-hat tension and why ``BoundarySourceSink`` dissolves it
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The earlier in-flight plan
(``.claude/plans/b52_boundary_residual_retype.md``) proposed typing the
matvec output boundary as :class:`BoundaryResidual`. That choice was
**rejected** for two reasons:

1. **It breaks consistency with the already-landed bulk.** The bulk
   ``.apply.bulk`` uses the source/sink leaf
   (:class:`AngularSourceSink`), **not** a residual, for operator
   outputs. Typing the boundary output as a residual would make the two
   halves of the same carve disagree on what an ``.apply`` output *is*.

2. **It creates a "two-hat" tension that the class gate cannot
   satisfy.** The realized boundary law
   :class:`~orpheus.sn.boundary_operator.SNBoundaryOperator` (:math:`B`)
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

   One operator cannot emit :class:`BoundaryResidual` for the matvec
   **and** :class:`BoundarySourceSink` for the SI rhs — the
   :class:`TimedFullField <orpheus.transport.timed_full_field.TimedFullField>`
   class gate (strict class identity:
   ``type(self.boundary) is not type(other.boundary)`` ⟹ ``TypeError``)
   throws on ``BoundaryResidual + BoundarySourceSink`` the moment the SI
   rhs tries to add :math:`B.\text{apply}` (a residual, under OPT-BR)
   to :math:`S.\text{apply}` and :math:`q_{\rm ext}` (sources). The
   variadic driver (:ref:`bc-extraction-variadic-driver`) makes this
   sharper than the retired fold: each gain's output is summed
   *individually*, so :math:`B`'s lone hat must be a source/sink for the
   rhs sum :eq:`bc-extraction-variadic-matvec` to close.

Choosing :class:`BoundarySourceSink` for **all** operator outputs
dissolves the two-hat: :math:`B` wears **one** hat (it always emits a
source/sink), and **both** sums close as homogeneous
:class:`BoundarySourceSink` sums —

.. math::
   :label: bc-extraction-two-hat-closed-sums

   \underbrace{(L+C).\text{apply} - S.\text{apply}
               - B.\text{apply}}_{\text{Krylov matvec}}
   \quad\text{and}\quad
   \underbrace{q_{\rm ext} + S.\text{apply}
               + B.\text{apply}}_{\text{SI rhs}}

both stay within the single :class:`BoundarySourceSink` class. This
needs **no SI-driver restructure** and **no partial-O.2**:
:class:`BoundaryResidual` stays reserved for the named
``from_balance`` composition exactly as
:class:`AngularResidual` waits on the bulk.

.. vv-status: bc-extraction-two-hat-closed-sums documented

A throwaway **decision instrument**
(``derivations/diagnostics/diag_b52_boundary_typing_decision.py``, the
B0 de-risk) proved on a 1-D reflective slab **and** a 2-D reflective
box that the OPT-BSS choice (``BoundarySourceSink`` for the matvec
output) closes both sums, while the OPT-BR choice
(:class:`BoundaryResidual` for the matvec output) throws the two-hat
``TypeError`` on the SI rhs.


Why the Krylov path is safe with a ``BoundarySourceSink`` matvec output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The matvec output never escapes scipy as a :class:`BoundarySourceSink`,
so the *solution* side stays :class:`BoundaryFlux`. The mechanism is the
flat round-trip:

* :meth:`TimedFullField.to_flat <orpheus.transport.timed_full_field.TimedFullField.to_flat>`
  ravels the composite to ``[bulk.values.ravel(), boundary.values]`` —
  a **type-agnostic** 1-D vector (the class of ``.boundary`` is erased).
* scipy's GMRES iterate is reconstructed via
  :meth:`TimedFullField.from_flat <orpheus.transport.timed_full_field.TimedFullField.from_flat>`,
  which rebuilds the boundary with
  ``replace(template.boundary, values=...)`` off the flux
  ``solution_template``. Because the template's boundary is a
  :class:`BoundaryFlux`, the reconstructed iterate's boundary is a
  :class:`BoundaryFlux`.

So the matvec's *internal* :class:`BoundarySourceSink` boundary class
lives only inside one ``op.apply`` call; the moment the result is
raveled and handed to scipy, the class is gone, and the iterate scipy
hands back is reconstructed as the flux type. The solve/iterate/trace
sites are therefore correctly **kept** :class:`BoundaryFlux`:
:meth:`CollisionOperator.solve <orpheus.sn.operator.CollisionOperator>`,
the boundary buffer of
:meth:`InvertibleOperator._solve_timed_full_field <orpheus.sn.operator.InvertibleOperator._solve_timed_full_field>`,
the cold-start ``initial_guess`` iterates
(``TimedFullField.zeros(..., boundary=BoundaryFlux, ...)``), the
converged traces, and the sweep's persistent boundary buffer.


The 13 retyped sites
~~~~~~~~~~~~~~~~~~~~

Thirteen sites (operator outputs + ``q_ext`` sources) flipped from
:class:`BoundaryFlux` to :class:`BoundarySourceSink`:

.. list-table:: B.5.2 retyped sites
   :header-rows: 1
   :widths: 34 38 28

   * - Module / symbol
     - Site
     - Emission
   * - :mod:`orpheus.sn.operator`
     - :meth:`_SpatialSweepDirection.apply <orpheus.sn.operator.StreamingOperator>`
     - bare-sweep boundary block
   * - :mod:`orpheus.sn.operator`
     - ``_compute_LpC`` (``m_boundary``)
     - :math:`L+C` boundary block
   * - :mod:`orpheus.sn.operator`
     - ``_compute_decomposition`` (``m_spat_boundary``, the
       **dual-emission twin** — flipped identically)
     - :math:`L+C` boundary block (decomposed)
   * - :mod:`orpheus.sn.operator`
     - ``_compute_decomposition`` (``M_angular`` zero)
     - auto-allocated zero
   * - :mod:`orpheus.sn.operator`
     - :meth:`StreamingOperator._apply_2d_cartesian <orpheus.sn.operator.StreamingOperator>`
     - 2-D boundary block
   * - :mod:`orpheus.sn.operator`
     - :meth:`CollisionOperator.apply <orpheus.sn.operator.CollisionOperator>`
     - bulk → bulk; boundary zero
   * - :mod:`orpheus.sn.scattering`
     - :meth:`ScatteringOperator.apply <orpheus.sn.scattering.ScatteringOperator>`
     - boundary zero
   * - :mod:`orpheus.sn.fission`
     - :meth:`FissionOperator.apply <orpheus.sn.fission.FissionOperator>`
     - boundary zero
   * - :mod:`orpheus.sn.boundary_operator`
     - :meth:`SNBoundaryOperator._apply_faces <orpheus.sn.boundary_operator.SNBoundaryOperator>`
       (``apply`` **and** ``apply_transpose``)
     - :math:`B\,\psi.\text{outflow}` on the inflow slots
   * - :mod:`orpheus.sn.solver`
     - ``_zero_within_group_fission`` (the ``F = ZeroOperator``
       codomain zero)
     - codomain zero
   * - :mod:`orpheus.sn.solver`
     - ``q_ext.boundary`` at the **3 source builds** (eigenvalue SI,
       eigenvalue Krylov, fixed-source SI / reconstruction)
     - prescribed inflow (zero for vacuum / reflective)

The change is **type-only**: :meth:`BoundarySourceSink.zeros_on <orpheus.transport.fields._bases.BoundaryField.zeros_on>`
and the per-face-view writes produce **bit-identical** ``.values`` —
only the wrapping role-type differs. The dead :class:`BoundaryFlux`
runtime imports were retired from the retyped sites.


.. _bc-extraction-operator-output-o2:

What remains for Wave O step O.2
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Wave O step **O.2a has landed the honest** :math:`L+C-S-B` **driver**
via the variadic couplings of :ref:`bc-extraction-variadic-driver`:
the transitional :math:`S + B` fold is **retired** and :math:`B` is now
a first-class coupling gain. What B.5.2 *still* leaves for the rest of
O.2 is the **residual column** of the role grid — and the adjoint
metric:

* :meth:`BoundaryResidual.from_balance <orpheus.transport.residuals.boundary_residual.BoundaryResidual.from_balance>`
  has **no operator-output consumer** yet. The honest variadic driver
  emits each gain's output as a :class:`BoundarySourceSink` (it never
  forms a typed field residual — the GMRES defect is the *flat*
  :math:`b - A\psi` on the raveled vector). The first consumer of
  ``BoundaryResidual.from_balance(lhs=ψ.inflow, rhs=B·ψ.outflow + q)``
  is the O.2 named-composition driver that types the affine boundary
  balance explicitly at the solver level.
* the :math:`|\Omega\cdot\hat n|\,w` :math:`G`-metric adjoint ``.H``
  (the boundary-weighted inner product for the transpose) — which is
  exactly why :math:`B` stays trace-typed as a separate gain
  (:ref:`bc-extraction-variadic-driver`), and
* **Gate-1.3** (the O.2 verification gate).

The direct-helper
:func:`_reflect_outflow_into_inflow <orpheus.sn.solver._reflect_outflow_into_inflow>`
also survives O.2a (the driver no longer routes through it, but the
final eigenvalue reconstruction sweep and the Gauss-Seidel variant
still do — :ref:`bc-extraction-two-routes`); the optional
:attr:`ScatteringOperator.domain` typing-completion tripwire remains a
documented seam (:ref:`bc-extraction-scope-future`).

Until the residual column is wired, :class:`BoundaryResidual` and
:class:`~orpheus.transport.residuals.angular_residual.AngularResidual`
are minted-but-consumerless role-grid completions — the correct status
the V&V audit reports for them.


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


.. _wavefront-flux-cochain:

The interior face-flux cochain — :class:`WavefrontFlux` (Wave O step #205 Phase 5)
==================================================================================

Wave O step #205 Phase 5 (`Issue #208
<https://github.com/deOliveira-R/ORPHEUS/issues/208>`_, commits
``478723d`` mint / ``992b0c0`` 2-D sweep / ``0e3e16c`` 2-D matvec,
2026-06-04) typed the SN wavefront sweep's **interior** cell-face
angular fluxes — historically the raw ephemeral numpy arrays ``psi_x``
``(N, ng, nx{+}1, ny)`` / ``psi_y`` ``(N, ng, nx, ny{+}1)`` — as a
named field :class:`~orpheus.transport.fields.wavefront_flux.WavefrontFlux`.
This is the **face-locus** sibling of the boundary-block typing of
:ref:`bc-extraction` (cell + trace) and the operator-output typing of
:ref:`bc-extraction-operator-output-typing`: where those typed the
*cell* state and the *operator outputs*, this types the *interior
face* state that the wavefront sweep propagates between them. It kills
the ``coding-elegance`` Pattern-3 anti-pattern (a flux-bearing tensor
with no type identity) and **names the trace operator the sweep
applies by hand**.


The native frame — discrete exterior calculus / cochains
--------------------------------------------------------

The native mathematical structure of the SN face fluxes (validated by
the cross-domain-attacker, frame memo
``field_role_typing_faceflux_frames.md``) is **discrete exterior
calculus**. The per-ordinate angular flux crossing faces is a primal
**1-cochain** — an assignment of a value to each oriented face:

.. (vv-status rationale) Structural framing of the SN face fluxes as a
   primal 1-cochain. This is a definitional / representational identity
   (the named-field typing), not a solver claim; the verifiable content
   is the bit-identity of the typed walk against the raw psi_x/psi_y
   walk, pinned by the octant-equivalence and Gate-K suites below, plus
   the 25 foundation tests of test_wavefront_flux.py.
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
     C^1_{\rm int} = \texttt{WavefrontFlux} & (\text{interior faces}), \\
     C^1_\partial  = \texttt{BoundaryFlux}  & (\text{domain-edge faces}),
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
       (:class:`~orpheus.transport.fields.boundary_flux.BoundaryFlux`)
     - the streaming sweep + sibling :math:`-B`
   * - **face**
       (:eq:`wavefront-cochain-biproduct`)
     - :math:`C^1_{\rm int}`
       (:class:`~orpheus.transport.fields.wavefront_flux.WavefrontFlux`)
     - :math:`C^1_\partial`
       (:class:`~orpheus.transport.fields.boundary_flux.BoundaryFlux`)
     - the trace operators :math:`\iota_*` / :math:`\iota^*`

:class:`BoundaryFlux` is the boundary summand at **both** loci — it is
literally the domain-edge faces, which are exactly the
boundary trace of the cell+trace state. The interior summand differs:
the cell biproduct carries the cell-centre flux, the face biproduct
carries the interior *face* flux. The boundary persists (it carries
the SI / Krylov iterate across calls); the interior is **ephemeral**
(rebuilt each sweep), so the two summands have different lifetimes —
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
write-back (the absorb). :class:`WavefrontFlux` names them:

.. list-table:: The typed trace operators
   :header-rows: 1
   :widths: 14 22 28 36

   * - Symbol
     - Method
     - Direction
     - Role in the biproduct
   * - :math:`\iota_*`
     - :meth:`~orpheus.transport.fields.wavefront_flux.WavefrontFlux.seed`
     - :math:`C^1_\partial.\text{inflow} \to C^1_{\rm int}`
       domain-edge slots
     - the **injection** :math:`\iota_\partial` of the biproduct
   * - :math:`\iota^*`
     - :meth:`~orpheus.transport.fields.wavefront_flux.WavefrontFlux.absorb`
     - :math:`C^1_{\rm int}` domain-edge slots
       :math:`\to C^1_\partial.\text{outflow}`
     - the **projection** :math:`\pi_\partial` of the biproduct

Both route through the single-source-of-truth
:meth:`_edge_slot <orpheus.transport.fields.wavefront_flux.WavefrontFlux>`
face-to-edge map, so the injection and projection cannot desync. A
third read accessor
:meth:`~orpheus.transport.fields.wavefront_flux.WavefrontFlux.edge_view`
exposes the domain-edge slot as a zero-copy view *without* copying into
a :class:`BoundaryFlux` — the 2-D matvec uses it to difference
``edge_view(face) − given`` when emitting its active-trace boundary
residual (so there is no hardcoded ``psi_x[:, :, 0, :]`` literal at the
call site).

The two biproduct laws follow, and are **provable** rather than
coincidental:

.. (vv-status rationale) The two biproduct identities — absorption =
   identity (project-after-inject) and projection-annihilates-the-
   strictly-interior (the off-diagonal block is zero). Structural laws
   of the biproduct; pinned by test_absorption_is_identity (slab /
   sphere / 2-D box) and test_pi_int_after_injection_is_zero_2d.
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
an *observation* under the raw-numpy seed/absorb; it is now a **named
biproduct law** pinned by ``test_absorption_is_identity`` across slab /
sphere / 2-D box. The second — :math:`\pi_{\rm int} \circ
\iota_\partial = 0` — is the biproduct's off-diagonal-zero condition:
injecting the boundary leaves the **strictly**-interior faces
(positions :math:`1 \ldots n{-}1` along each axis) untouched at zero
(``test_pi_int_after_injection_is_zero_2d``).


Why the interior cochain is flux-only (single role)
---------------------------------------------------

:class:`WavefrontFlux` carries the **flux** role only — unlike the
boundary trace, it has no source / residual leaves. The reason is
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
strictly-interior faces carry no such balance, so :class:`WavefrontFlux`
is flux-only by construction — ``illegal-states-unrepresentable``: there
is no ``InteriorFaceResidual`` to mistype an interior face as.


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
     - :class:`~orpheus.transport.fields.wavefront_flux.WavefrontFlux`
       (**#205 Phase 5**)
     - — (off-diagonal of :math:`L_{\rm oct}`, no source role)
     - — (no cell-balance defect on interior faces)
   * - **face — boundary**
       (1-cochain :math:`C^1_\partial`)
     - persistent (carries the iterate trace)
     - :class:`~orpheus.transport.fields.boundary_flux.BoundaryFlux`
     - :class:`~orpheus.transport.source_sinks.boundary_source_sink.BoundarySourceSink`
       (B.5.2)
     - :class:`~orpheus.transport.residuals.boundary_residual.BoundaryResidual`
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

:class:`WavefrontFlux` stores a **single flat backing buffer**
(``space.layout.total_size``); the per-axis face fields are zero-copy
reshape views (:meth:`~orpheus.transport.fields.wavefront_flux.WavefrontFlux.face`).
The cross-domain-attacker **rejected** a per-face Python object on three
independent grounds:

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

This mirrors :class:`BoundaryFlux`, which already uses
``flat-buffer + FaceLayout + face_view``. The interior space
:class:`~orpheus.numerics.spaces.interior_face_space.InteriorFaceSpace`
is the **layout-on-space** (A.5) sibling of
:class:`~orpheus.numerics.spaces.trace_space.TraceSpace` minus the
``omega_dot_n`` directional selectors (the interior cochain has no
inflow/outflow partition — it is flux-only). It is axis-parametric: the
:meth:`~orpheus.numerics.spaces.interior_face_space.InteriorFaceSpace.interior_layout`
builder takes the axis count as a parameter (one face-normal field per
active axis), so a future 3-D Cartesian wavefront sweep is a
*parameter* (``axes=(0,1,2)``), not a new code path — the 3-D-readiness
the cochain frame's dimension-agnostic structure licenses
(``feedback_unify_after_two_instances``: 1-D + 2-D are two working
instances, 3-D is the validating third).


.. _wavefront-flux-honest-scope:

Honest scope — a representation win, NOT a speed/rate/parallelism win
---------------------------------------------------------------------

.. warning::

   :class:`WavefrontFlux` is a **representation / elegance win only**.
   It does **NOT**:

   * change the asymptotic cost of the sweep,
   * recover the source-iteration convergence rate,
   * enable parallelism on its own.

   The cross-domain-attacker stated this as an *absence*, not a hedge:
   the wavefront DAG already gives the optimal sweep schedule, and the
   typing relocates no arithmetic. The seed/absorb stays an **inherent
   cheap** :math:`O(\text{boundary faces})` copy — negligible against
   the :math:`O(\text{volume})` sweep — at the **persistent-boundary /
   ephemeral-interior lifetime split**. True zero-copy between the two
   summands is *precluded* (and unnecessary) because
   :class:`BoundaryFlux` persists across SI iterations while
   :class:`WavefrontFlux` is rebuilt each sweep.

The payoff is the **type**: a named field, typed :math:`\iota_*` /
:math:`\iota^*`, code that reads like the cochain math, and
illegal-states-unrepresentable (the flux-only constraint of the
interior locus is enforced by there being no interior-face residual
leaf). It is also the **clean substrate** the SI Gauss–Seidel recovery
lands on: with typed :class:`WavefrontFlux` + :class:`BoundaryFlux` +
:math:`\iota^*`, the ``(octant × face)`` reflective-graph G-S schedule
becomes an explicit typed composition (``sweep octant → ι* absorb → −B
reflect → ι_* seed next octant``) rather than an implicit buffer-timing
dance — that recovery, on this substrate, is where the actual
convergence-rate win is sought (a separate, research-tagged step).


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
   * - **25 foundation tests**
       (``tests/transport/fields/test_wavefront_flux.py``)
     - units / class identity / field+views / the two biproduct laws /
       the round-trip pin + L11 negative control / axis-parametricity
       (1-D / 2-D / 3-D one path)
     - all ``@pytest.mark.foundation`` (software invariants, no theory
       ``:label:``)
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

   :class:`WavefrontFlux` types the **2-D** wavefront sweep + matvec
   only. The **1-D** sweep is a parallel-prefix **scan** (a different
   fold): its interior fluxes are transient chain-ordered scan output
   (``(nx, K, ng)``, no persistent interior buffer at all), the
   structural antithesis of the cochain ``(N, ng, nx{+}1)``. Forcing
   :class:`WavefrontFlux` into the 1-D scan would be a wrong-fit (it
   risks the L16 ``cumprod`` efficiency and multiplies concepts). The
   *type* is built axis-parametric (it accepts ``axes=(0,)`` through
   the same code path — the realised 1-D benefit), but the 1-D
   *implementation* keeps its scan fold. The one shared seam (the
   boundary-trace exchange + DD-closure averaging) unifies cleanly only
   when a future ``nd_foundation`` session re-expresses **both** folds
   as one :math:`d`-generic walk, under the hard constraint that the
   collapse must not regress the :math:`d=1` scan's parallel-prefix
   efficiency. A future 3-D *Cartesian* sweep IS a wavefront and uses
   :class:`WavefrontFlux` directly (``axes=(0,1,2)``).


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
:class:`~orpheus.numerics.spaces.trace_space.InflowTraceSpace` and
:class:`~orpheus.numerics.spaces.trace_space.OutflowTraceSpace`, which
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
:meth:`InflowTraceSpace.from_mesh_and_quadrature`, which builds
the mask from the spatial mesh's face-normal table and the
quadrature's direction-cosine arrays. Every :class:`Mesh1D` coord
system (``CARTESIAN`` / ``SPHERICAL`` / ``CYLINDRICAL``) shares the
same ``("left", "right")`` face structure with the radial axis as
the outward normal — :class:`GaussLegendre1D` is the shared
quadrature, with ``mu_x`` as the direction cosine along that axis.
The :class:`Mesh2D` factory supports ``coord=CARTESIAN`` only;
2-D cylindrical (axisymmetric :math:`(r, z)`) raises
:class:`NotImplementedError` and will until a 2-D cylindrical SN
sweep ships in ORPHEUS (none exists today).

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
       :math:`|\Omega\cdot\hat n|\,w` :math:`G`-metric (codomain inner
       product of :math:`L`'s boundary-trace block) lives HERE and is
       reused by every posing.
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
       BiCGSTAB. Diffusion: FD solve. Inverts *whatever*
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
**no** :math:`(L, S, F)` split at all — its
:meth:`solve_fixed_source <orpheus.numerics.eigenvalue.EigenvalueSolver.solve_fixed_source>`
is one BiCGSTAB on a *monolithic* collision-probability matrix; the
factor :math:`(L-S)^{-1}` does not exist as a separable object.
Splitting the posing into

* **2a — role assignment + μ-map** (pure data: a posing-table row), and
* **2b — loss-operator realization** (the method's concrete assembly),

lets :math:`2a \circ 2b \circ 3 \circ 4` compose cleanly across every
family. The key consequence:
:class:`~orpheus.numerics.iteration.KEigenvalue(L, S, F)` is the
operator-triple **2b realization** — NOT a problem-type layer. Treating
the operator triple as a "problem type" was the conflation the
bifurcation removes.

**The variadic driver IS the posing/resolvent boundary made explicit.**
The Layer-3 SN resolvent
(:class:`~orpheus.numerics.iteration.SourceIteration` /
:class:`~orpheus.numerics.iteration.KrylovAcceleration`) is now
**problem-type-agnostic**: it consumes
:math:`\text{Driver}(L_{\rm resolvent},\,*\text{gains})` and never asks
which gain plays which posing role. *Which* leaves are gains is the
2a decision — for the SN k-row the gains are exactly the
:math:`A_{\rm loss}` couplings :math:`S` and :math:`B` (fission
:math:`F` is the eigen-operator :math:`M`, not a within-group gain; it
enters as :math:`q_{\rm ext}`). The retired fixed :math:`(L, S, F)`
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
:eq:`tensor-product-adjoint-distributivity` and the
:math:`|\Omega\cdot\hat n|\,w` :math:`G`-metric ``.H`` of
:ref:`bc-extraction-operator-output-o2`). Because forward and adjoint
share the spectrum, :math:`\lambda` — and therefore the
:math:`\mu \to` physical map — is **unchanged**; only the operators are
daggered. The adjoint slots in at 2a with **zero new engine
machinery**: it is a row, not a layer. (The first-draft architecture's
instinct to make adjoint a separate "mode" is the same conflation the
2a/2b split removes.)

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
  **early**, building it as :math:`(L-S)^{-1}` from the operator triple
  via an inner :class:`~orpheus.numerics.iteration.SourceIteration`.

The late-bound layer is **strictly more general**: it admits *both* the
operator-triple resolvent (SN, MoC) *and* the **monolithic-matrix
resolvent** (CP, diffusion, homogeneous) that has no :math:`(L, S, F)`
factorization. The early-bound layer can only express methods whose
resolvent factors as :math:`(L-S)^{-1}` from a triple — strictly
narrower. The general layer cannot be expressed in terms of the narrow
one without *manufacturing fictitious* :math:`L`, :math:`S` operators
for CP/diffusion (which have no sweep). Therefore the **Protocol layer
is canonical** and the **triple layer is a specialization that adapts
into it**.

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
     - **early** — :math:`(L-S)^{-1}` from the triple
   * - Admits
     - SN, MoC, CP, diffusion, homogeneous (any Protocol implementer)
     - SN / MoC only (needs an :math:`(L,S,F)` triple)
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
:math:`= (L-S)^{-1} q` via the warm-started inner
:class:`~orpheus.numerics.iteration.SourceIteration`, the injected
:math:`k`- and production-estimators, and the :math:`\ge 3`-iteration
:math:`dk`/:math:`d\phi` convergence test — then
:meth:`solve <orpheus.numerics.iteration.KEigenvalue.solve>` simply
calls ``power_iteration(self, max_iter=self.max_outer)``. SN production
(:func:`~orpheus.sn.solver.solve_sn`), CP, diffusion, MoC, and
homogeneous all drive the same loop directly via the Protocol; the
``KEigenvalue`` adapter is for callers who *have* a natural
:math:`(L,S,F)` triple and want to skip writing a full solver class.

.. note::

   The "``KEigenvalue`` regresses :math:`P_\ell` (anisotropic
   scattering)" objection from the migration era dissolves under this
   framing: :class:`~orpheus.numerics.iteration.SourceIteration` is
   type-agnostic and **angular-capable** — it routes the RHS through
   the ravellable protocol, so a typed
   :class:`~orpheus.sn.scattering.ScatteringOperator` acting on a
   :class:`~orpheus.transport.timed_full_field.TimedFullField` carries
   :math:`P_\ell` correctly. The observed regression was a property of
   an L1 *test adapter* that collapsed angular flux to scalar between
   outer iterations (dropping the angular moments), not of
   ``KEigenvalue``. The decisive — and sufficient — reason
   ``KEigenvalue`` cannot be the universal engine is the
   CP/diffusion/homogeneous **no-triple** fact alone.


The metric lives at the leaf, not the posing
---------------------------------------------

The boundary :math:`|\Omega\cdot\hat n|\,w` :math:`G`-metric is the
*codomain inner product* of :math:`L`'s boundary-trace block (the Wave
O dagger-biproduct, :ref:`bc-extraction-operator-output-o2`). It is
**intrinsic to the streaming leaf** — the leaf carries its ``domain`` /
``codomain`` :class:`~orpheus.numerics.space.FunctionSpace`\ s with
their ``inner_product_weights``, and the ``.H`` adjoint reads them. It
is **NOT a posing-layer concern**: posing *arranges* leaves; the leaves
already *know* their metric. Consequently the adjoint posing row gets
the correct :math:`G`-weighted transpose for free — the dagger functor
applied to leaves that already carry the metric — which is precisely
why the adjoint row adds no new machinery.


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

