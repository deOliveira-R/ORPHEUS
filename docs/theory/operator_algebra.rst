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

