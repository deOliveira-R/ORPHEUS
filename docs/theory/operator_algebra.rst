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


.. _bc-extraction-two-routes:

The two :math:`-B` delivery routes (and their O.2 collapse)
-----------------------------------------------------------

The same :math:`-B` coupling reaches the iteration two ways, both
calling the **identical**
:class:`~orpheus.sn.boundary_operator.SNBoundaryOperator` (single
source of truth, Cardinal Rule 2):

.. list-table:: The two delivery routes for :math:`-B`
   :header-rows: 1
   :widths: 22 40 38

   * - Route
     - Mechanism
     - Used by
   * - **Driver fold**
       (:meth:`SNSolver._scattering_with_boundary_op <orpheus.sn.solver.SNSolver._scattering_with_boundary_op>`)
     - Folds :math:`S + B` into the **subtracted** :math:`S` argument
       of the iteration drivers: the matvec
       :math:`(L+C).\text{apply} - (S+B).\text{apply} - F.\text{apply}
       = (L+C-S-F-B)`. :math:`B` cannot join the :math:`L+C`
       preconditioner because
       :class:`~orpheus.numerics.operator.OperatorSum` does **not**
       propagate ``CAP_SOLVE`` — :math:`L + C - B` would strip the
       ``.solve`` (sweep) the Krylov preconditioner / SI splitting
       needs.
     - The eigenvalue SI inner driver
       (:meth:`SNSolver._solve_source_iteration <orpheus.sn.solver.SNSolver._solve_source_iteration>`),
       the eigenvalue Krylov inner
       (:meth:`SNSolver._solve_krylov <orpheus.sn.solver.SNSolver._solve_krylov>`),
       and the fixed-source Krylov
       (:func:`_solve_fixed_source_krylov <orpheus.sn.solver._solve_fixed_source_krylov>`).
       :math:`B\,\psi.\text{outflow}`
       rides in :math:`\text{rhs.boundary}`, which the bare
       ``(L+C).solve`` sweep reads as the inflow seed.
   * - **Direct helper**
       (:func:`_reflect_outflow_into_inflow <orpheus.sn.solver._reflect_outflow_into_inflow>`)
     - Fills each face's inflow slots with
       :math:`B\,\psi.\text{outflow}` in place on the boundary buffer,
       via the same :class:`SNBoundaryOperator`, before the bare
       sweep.
     - The DIRECT loops that have no driver to fold into: the direct
       fixed-source SI
       (:func:`_solve_fixed_source_si <orpheus.sn.solver._solve_fixed_source_si>`)
       and the final eigenvalue reconstruction sweep in
       :func:`solve_sn <orpheus.sn.solver.solve_sn>`. Guarded to
       1-D via ``sn_mesh.reduced is not None`` (see scope below).

The two routes **collapse at Wave O step O.2**: when the iteration
drivers take the whole loss operator :math:`L+C-S-F-B` directly (on
the direct-sum carrier), the direct loops route through the driver and
:func:`_reflect_outflow_into_inflow` retires (its documented removal
trigger).

**The O.2 forcing function.** The :math:`S + B` fold type-checks
**only because** :attr:`ScatteringOperator.domain` is ``None`` (it
predates function-space tagging). The
:class:`~orpheus.numerics.operator.OperatorSum` domain-compatibility
check fires only when both operands declare non-``None`` domains that
differ; with :math:`S` untagged it skips. The moment the typing wave
gives :math:`S` a (bulk-space) domain, the fold throws
``IncompatibleOperatorComposition`` at construction — :math:`S` (bulk
space) and :math:`B` (trace space) live on different function spaces.
That throw **is** the O.2 signal to land the honest composition: the
drivers consuming the whole :math:`L+C-S-F-B` loss operator on the
direct-sum carrier, retiring the fold and the helper together. The
fold is documented in code as an O.2 *tripwire*, not a permanent
design.


.. _bc-extraction-scope:

Scope — what is bare, what is not (O.4b)
----------------------------------------

O.4a.2 made the **1-D** sweep bare (slab / sphere / cylinder). The
**2-D Cartesian wavefront sweep is NOT yet bare** — it still applies
``bc`` internally (the :ref:`sweep-octant-dependency-graph` L7-trap
fix applies the BC once per octant per axis, but it is still applied
*inside* the sweep). That extraction is deferred to **O.4b**;
:class:`~orpheus.sn.boundary_operator.SNBoundaryOperator` is not yet
wired for the 2-D trace.

The dispatch is guarded by a **single predicate** so the two paths
cannot drift: ``sn_mesh.reduced is not None`` is the **same** predicate
:func:`~orpheus.sn.sweep.transport_sweep` uses to select its
bare-vs-bc-in-sweep body, and the **same** predicate the direct-helper
guards
(:func:`_solve_fixed_source_si <orpheus.sn.solver._solve_fixed_source_si>`,
:func:`solve_sn <orpheus.sn.solver.solve_sn>`) check before calling
:func:`_reflect_outflow_into_inflow <orpheus.sn.solver._reflect_outflow_into_inflow>`.
1-D goes through the bare sweep + :math:`-B`; 2-D stays on the
unchanged bc-in-sweep wavefront until O.4b.


.. _bc-extraction-numerical-evidence:

Numerical evidence
------------------

The extraction is verified by three independent grounds (per the
``vv-principles`` skill's three pillars and the bit-identity
vs principled-equivalence gate).

**1. Vacuum bit-identity.** With :math:`B = 0` the fold
:math:`S + B \equiv S` exactly, so the vacuum path is **bit-identical**
to the pre-extraction matvec. Verified by:

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

