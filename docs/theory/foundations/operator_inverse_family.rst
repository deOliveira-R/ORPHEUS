.. _operator-inverse-family:

===========================
The Operator Inverse Family
===========================

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
      concept: operator_inverse_family
      role: "how the operator algebra's inverse A⁻¹ is realized and materialized (#226 taxonomy)"
      covers: [driver-applied sweep inverse, Green preconditioned-splitting inverse, dense materialising inverse, sparse assembly functor]
      depends_on: [operator_algebra]
      status: "extracted from operator_algebra.rst; content verbatim, provisional header"


This page is the **inverse family** of the operator algebra: the four
realizations of "apply :math:`A^{-1}` / materialize :math:`A`" that the #226
operator taxonomy separates. The operators being inverted — the loss
composite :math:`A = L + C - S - N_{2n} - B`, and its invertible
sub-composite
:math:`L + C` whose inverse **is** the transport :term:`sweep` — are developed in
:doc:`/theory/foundations/operator_algebra`. This page picks up at the
inverse itself and follows it through four distinct realizations:

- the **solver builds / driver applies** split — the posing layer builds
  :math:`A^{-1}` once and the iteration primitive *applies* it
  (:ref:`inverse-application-driver`);
- the **Green operator**, the inverse of a named split :math:`M - N`,
  minted by the Strategy with the realisation it chose, and the first
  *iterative* member of the family (:ref:`green-operator`);
- the **dense materialising inverse** — ``as_matrix`` out of the operator
  category, factored once into a
  :class:`~orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator`
  (:ref:`matrix-inverse-operator`);
- the **assembly functor** — structural sparse emission, a homomorphism
  out of the operator algebra living downstream of it, and the second
  algorithm for the matrix ``as_matrix`` probes
  (:ref:`operator-algebra-assembly-axis`).


.. _inverse-application-driver:

The solver builds the inverse; the driver applies it
====================================================

Steps 1–2 reified the inverse *operators* — the schedule-triangular
:class:`~orpheus.sn.operators.sweep_operator.SweepOperator` (step 1, the
WDD sweep as :math:`(L+C)^{-1}`), the reified splitting matrix
:math:`M=(L+C)-B_{\rm lower}` (step 2, :ref:`si-gauss-seidel-reification`
in :doc:`/theory/methods/sn/cartesian_multid`), and the windowed composition
:math:`P\circ A^{-1}` (step 2, :ref:`windowing-retyped` in
:doc:`/theory/methods/sn/cartesian_multid`). **Step 3
retypes the driver–inverse boundary.** The solver (posing) layer builds
the inverse *once*, and the iteration primitive
(:class:`~orpheus.numerics.iteration.SourceIteration`) *applies* it:

.. (vv-status rationale) The SI driver update
   ψ_{n+1} = A_inv.apply(q_ext + Σ_i g_i ψ_n, initial_guess=ψ_n), with
   A_inv = A⁻¹ built by the solver (#226 step 3). Governing iteration
   (structural); the SI fixed point is exercised by the source-iteration
   solver suites.
.. vv-status: inverse-driver-si-update documented

.. math::
   :label: inverse-driver-si-update

   \psi_{n+1} \;=\; A_{\rm inv}.\mathrm{apply}
   \Bigl(q_{\rm ext} + \textstyle\sum_i g_i\,\psi_n,\;
   \ \mathrm{initial\_guess}=\psi_n\Bigr),
   \qquad
   A_{\rm inv} \;=\; A^{-1}\ \text{(built by the solver)} .

Concretely, the S\ :sub:`N` source iteration's ``A_inv`` is the
preconditioner of the Strategy's Green: the
:class:`~orpheus.numerics.green_operator.GreenOperator` builds
:math:`M^{-1}` once, its
:class:`~orpheus.numerics.green_operator.Richardson` realisation derives
the step operator from it (the sweep itself, or the windowed product
``P @ M.inverse()`` on 2-D Cartesian), and the Green's relaxation face
hands that step operator to the ``SourceIteration`` constructor as its
**first argument** (:ref:`green-strategy-mints`). This closes the #226
steps-1–3 arc: the
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

.. (vv-status rationale) The canonical SupportsSeededApply signature
   apply(rhs, *, initial_guess=None) → codomain. A notation / type-contract
   definition (a static annotation target, deliberately not
   runtime_checkable).
.. vv-status: seeded-apply-signature documented

.. math::
   :label: seeded-apply-signature

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
       starting-direction flux from the previous iterate's :term:`per-ordinate <ordinate>`
       :math:`\psi`; ``apply`` delegates to ``inner.solve(rhs,
       initial_guess=...)`` which reads it.
   * - :class:`~orpheus.sn.operators.sweep_operator.SweepOperator`
       (2-D Cartesian DD)
     - accepted, ignored
     - The :term:`diamond-difference <diamond difference>` wavefront is a direct forward substitution
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
**unconditionally** (no per-call decision) is what let the driver be
shape-agnostic — pinned at the time by the always-on Mode-11 spy
``test_seed_threading_spy.py`` (call *k*'s ``initial_guess`` must equal
call *(k−1)*'s return, by value; **retired in #280 2.5c** together with
the vestigial ``initial_guess`` threading it guarded). The live pins on
the seed contract are the strict no-forwarding spy in
``tests/gates/sn/operators/test_inverse_operator_equivalence.py`` and the
curvilinear value catcher in ``tests/gates/sn/eigenvalue/test_keff_curvilinear.py``.


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


The narrowing boundary and the structural resolution
----------------------------------------------------

``SupportsSeededApply`` and ``SupportsInverse`` are not
``runtime_checkable`` (an ``isinstance`` against a Protocol carrying only
a method member is a false-positive machine). The narrow from "an
invertible operator" to "a seeded-apply operator" therefore lives in
**one** place,
:func:`~orpheus.numerics.iteration.seeded_inverse` (spelled
``_seeded_inverse`` until #276 A4 promoted it — a third consumer, the
daggered adjoint :class:`~orpheus.numerics.iteration.SourceIteration`,
had joined the two below and cross-module callers were already importing
the underscore name):

.. code-block:: python

   def seeded_inverse(A):                        # numerics/iteration.py
       if not invertible(A):                     # the guard lives HERE
           raise NotInvertible(...)
       inv = A.inverse()
       if accepts_seed(inv):                     # TypeGuard: the ONE capability
           return inv                            # canonical seeded apply
       return _SeededExactApply(inv)             # accept-and-drop wrapper

Two things about that body are load-bearing. First, **the runtime
precondition is internal, not the caller's**: ``seeded_inverse`` runs its
own :func:`~orpheus.numerics.operator.invertible` check and raises
:class:`~orpheus.numerics.operator.NotInvertible`, so a caller's
``is_invertible`` guard is a *domain-message* courtesy, never the thing
that makes the call sound. Second, **there is no** ``cast``: both
narrowings are checked bridges (Design C) — ``invertible`` is a
``TypeGuard`` whose runtime predicate *is* the static permission for the
``.inverse()`` call, and
:func:`~orpheus.numerics.iteration.accepts_seed` is a ``TypeGuard`` that
reads the inverse's structure, never its signature: a wrap-delegate
member (``isinstance`` on
:class:`~orpheus.numerics.operator.InverseWrapMixin`), or a one-term
combination :math:`c\,w` whose word accepts a seed. The seedless branch
wraps in ``_SeededExactApply``, which accepts and drops the seed. That
branch is reached only by an exact inverse whose ``apply`` has no seed
slot (the identity, a permutation, a scaled identity), and only there is
"an exact inverse has nothing to seed" the whole truth; a scaled inverse
over a seeded word keeps its own seed law (:ref:`green-seed-law`).

Since package S4 of the numerics untangling the one production caller is
the :class:`~orpheus.numerics.green_operator.GreenOperator` constructor,
which builds :math:`M^{-1}` once as ``implicit_inverse``; every driver
that applies an inverse is built through a Green, so the narrowing (and
its guard) is single-sourced (`[M]` 2026-10-21, ``git grep -n
"seeded_inverse(" -- orpheus`` in the S4 worktree: 1 call outside
docstrings, ``green_operator.py``).

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
entirely. The narrowing ``seeded_inverse`` performs is therefore now a
checked bridge over an *already-guaranteed* shape rather than a
hoped-for one, and every wrap-delegate ``.inverse()`` a posing layer
reaches through it carries the seeded signature by construction.

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
rebuilt fresh objects). The factors stay reachable as ``.inner.factors``
(the flattened word, since package S2; the binary ``.inner.a`` /
``.inner.b`` retired with it).

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
     - :class:`~orpheus.numerics.elementary_operators.PermutationOperator` (inverse *is*
       a permutation), :class:`~orpheus.numerics.operator.IdentityOperator`
       (self-inverse), the one-term
       :class:`~orpheus.numerics.operator.OperatorSum` :math:`\alpha L`
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

* :class:`~orpheus.numerics.iteration.KEigenvalue` takes the Strategy's
  Green of the loss split, ``KEigenvalue(posing, green, …)``, and builds
  its inner once through the Green's relaxation face,
  ``green.driver(max_iter=max_inner, tol=inner_tol)``
  (:ref:`green-three-faces`). The invertibility of :math:`A` is the
  Green's construction obligation: a non-invertible implicit part raises
  :class:`~orpheus.numerics.operator.NotInvertible` when the Green is
  built, before ``KEigenvalue`` exists.
* :class:`~orpheus.numerics.iteration.KrylovAcceleration` keeps the
  **forward** :math:`A` (its GMRES matvec is
  :math:`A\cdot - \sum_i g_i\cdot`) and takes its left preconditioner
  :math:`P \approx \bigl(A - \sum_i g_i\bigr)^{-1}` as an operator the
  caller states, keyword-only with no default; plain GMRES is
  ``preconditioner=IdentityOperator()``.  The SN within-group builder
  states the sweep, ``(L+C).inverse()``, the transport-corrected
  preconditioner of Adams & Larsen 2002 §III, or the DSA-corrected
  :math:`(I + \mathcal{C})(L+C)^{-1}`; why there is no default is
  recorded at :ref:`sn-krylov-sweep-preconditioner`.  (At step 3 the
  driver's default-preconditioner fallback was rewired from the old
  ``CAP_SOLVE`` probe, with a ``# type: ignore``, to the honest
  ``A.is_invertible`` test plus ``seeded_inverse(A).apply``; the default
  itself retired on 2026-10-05, ``c1dc0c3d``.)


Verification — the seed spy and the windowed×G-S corner
-------------------------------------------------------

The rewire was pinned by three structural gates (all ``-O``-proof —
``pytest.fail`` / ``np.testing.assert_*``, never a bare ``assert``;
``test_seed_threading_spy.py`` was **retired in #280 2.5c** with the
vestigial ``initial_guess`` threading it guarded — the live seed-contract
pins are ``tests/gates/sn/operators/test_inverse_operator_equivalence.py`` and
``tests/gates/sn/eigenvalue/test_keff_curvilinear.py``):

.. list-table:: Step-3 regression gates
   :header-rows: 1
   :widths: 42 58

   * - Gate
     - What it pins
   * - ``test_seed_threading_spy.py``
       (foundation / sentinel; vv Mode 11)
     - Every inner solve's ``initial_guess`` equals the previous iterate's
       return, by value. **Route-invariant across the rewire**: it wraps
       :meth:`StreamingCollisionOperator.solve <orpheus.sn.operators.streaming.StreamingCollisionOperator.solve>`
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
       :class:`~orpheus.sn.operators.streaming.StreamingCollisionOperator` forward
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

The Green operator — the inverse of a named split
==================================================

Steps 1–3 delivered the *exact/direct* inverse family: the leaf
:class:`~orpheus.numerics.operator.InverseOperator` (a pointwise division)
and the schedule-triangular
:class:`~orpheus.sn.operators.sweep_operator.SweepOperator` (the WDD sweep
as :math:`(L+C)^{-1}`), both invertible in a single pass. **Step 4 adds the
first *iterative* member** — the
:class:`~orpheus.numerics.green_operator.GreenOperator`, built from a named
split :math:`A = M - N` as ``GreenOperator(M, N)``. It is the discrete
Green's function of the split operator: the object whose columns are the
flux responses to unit point sources.

.. math::
   :label: green-split-inverse

   \texttt{GreenOperator}(M, N) \;=\; (M - N)^{-1},
   \qquad M \text{ invertible},
   \qquad x \;\leftarrow\; M^{-1}\bigl(q + N x\bigr) .

A sum has no inverse of its own (:eq:`combination-inverse`): the
combination ``(L + C) - S`` reports ``is_invertible = False`` in every
spelling order, and its ``inverse()`` raises
:class:`~orpheus.numerics.operator.NotInvertible` naming
``GreenOperator(M, N)``. The inverse of a sum exists only once someone has
said which part is implicit, and that choice is the Strategy's
(``Splitting.from_schedule`` labels the loss terms), never the spelling's
(ruling 50 of ``.claude/plans/numerics_untangling.md``). ``inner`` is the
combination ``M - N`` itself, so ``G.inverse()`` returns it, and it equals
the sum spelled apart: ``GreenOperator(lc, S).inverse() == (lc - S)``.
:math:`M` must be invertible, and the constructor refuses a non-invertible
:math:`M` with ``NotInvertible`` naming its type.

.. implements:: green-split-inverse
   :by: orpheus.numerics.green_operator.GreenOperator

   **Implemented by** the class: the constructor's guard on :math:`M`, the
   forward ``inner = M - N``, :math:`M`'s own seeded inverse built once
   (``implicit_inverse``), from which a Richardson realisation derives its
   preconditioner, and the driven gains derived from :math:`N`'s terms
   (:ref:`green-realisations`).


The convergent :math:`A`-preconditioned splitting
-------------------------------------------------

A general operator sum has **no operand-wise inverse**: :math:`(A+B)^{-1}`
is not a function of :math:`A^{-1}` and :math:`B^{-1}`
(Sherman–Morrison–Woodbury applies only under low-rank structure — absent
from the bulk blocks :math:`L`, :math:`C`, but exactly what the boundary
block :math:`B` has: :ref:`the scoped statement <smw-low-rank-exception>`,
Issue #300), and —
unlike the schedule-triangular ``(L+C)`` — there is no substitution order
that solves it in one pass. What an operator *does* have, once it is
**split** as :math:`A - B` with an invertible implicit part :math:`A` (the
:math:`M` of :eq:`green-split-inverse`, the gains carried with their
signs in :math:`B`, the :math:`N`), is a convergent iteration.
Precondition by :math:`A^{-1}`:

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
:math:`k`-times-rescattered contribution (:cite:`LewisMiller1984` §3.2;
:cite:`AdamsLarsen2002` §II). Its partial sums are **exactly** the Richardson /
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

.. (vv-status rationale) The scattering-ratio contraction bound
   ρ((L+C)⁻¹S) ≤ max Σ_s/Σ_t = c < 1 (Adams & Larsen 2002 §II).
   Literature-transcribed mathematical identity, matching the sentineled
   green-neumann-series / green-splitting-iteration siblings.
.. vv-status: green-scattering-ratio-bound documented

.. math::
   :label: green-scattering-ratio-bound

   \rho\bigl((L+C)^{-1}S\bigr)
   \;\le\; \max_{\rm cell}\ \frac{\Sigma_s}{\Sigma_t} \;=\; c \;<\; 1 ,

guaranteed below unity for any absorbing (physical) medium
(:cite:`AdamsLarsen2002` §II: :math:`\rho = c`). The convergence is thus a
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
   G^{\mathsf T}\mathcal G` carrying the angular Gram metric. That
   adjoint exists on the object: ``G.dual()`` is the Green of the dual
   split, ``GreenOperator(M.dual(), N.dual())``, ``G.H`` is the base's
   ``♯ ∘ G.dual() ∘ ♭``, and by the law
   :math:`((M - N)^{-1})^{*} = (M^{*} - N^{*})^{-1}` it equals the Green of
   the adjoint split, which the honesty gate builds as its independent
   witness (:ref:`operator-adjoint-two-questions`). Both need the implicit
   part's inverse to transpose, which on S\ :sub:`N` is the 1-D reverse scan; on
   a multi-dimensional mesh the reverse scan is not realised
   (`#616 <https://github.com/deOliveira-R/ORPHEUS/issues/616>`_; the deferral ledger of
   :ref:`loss-rep-orientation-two-frames`), so there
   ``G.is_adjointable`` reads ``False`` and both refuse. This
   reciprocity gate does not use them: it builds :math:`G^{\mathsf T}`
   by hand from the transposed operands, so it stays independent of the
   object it would otherwise test.

The family is keyed by **which mathematical object the inverse is, not by
the algorithm that realizes it**. A Richardson-realized Green and a
GMRES-realized Green are the *same* Green operator — which is why
``KrylovInverseOperator`` was **rejected** as a sibling name: "Krylov"
names the orthogonal *realization* axis, not the object.


The split is named, never read off a spelling
---------------------------------------------

Which part of a loss is implicit decides the algorithm, so it is a choice
someone makes, and the type records it. The edges, as gated on the
heterogeneous 2-group vacuum slab of
``tests/gates/sn/operators/test_green_operator_sn.py``:

.. list-table:: How a split routes its inverse
   :header-rows: 1
   :widths: 26 26 48

   * - Spelling
     - Inverse
     - Behaviour
   * - ``build_streaming_collision(problem, mat_xs)``, the atom
       :math:`L + C`
     - ``lc.inverse()``, a
       :class:`~orpheus.sn.operators.sweep_operator.SweepOperator`
     - The direct sweep, carried by the atom.
   * - ``L + C`` or ``C + L``, the plain combination (equal)
     - none
     - ``is_invertible`` is ``False``; ``inverse()`` raises
       ``NotInvertible`` naming ``GreenOperator(M, N)``.
   * - ``A_loss = lc - S``, and ``(-S) + lc`` (the same combination)
     - none of its own; ``GreenOperator(lc, S)``
     - The physical splitting: preconditioner the sweep, gain :math:`S`.
       Converges (:math:`c < 1`), and ``GreenOperator(lc, S).inverse() ==
       A_loss``.
   * - ``GreenOperator(C, -L)``
     - a Green that constructs
     - A legal split whose implicit part happens to be invertible — the
       algebra cannot read :math:`\rho(C^{-1}L) > 1`. Its
       collision-preconditioned Richardson **diverges**, and ``apply``
       raises :class:`~orpheus.numerics.green_operator.ConvergenceFailure`
       **loudly**. Same operator as the physical split, a different
       algorithm by the choice of :math:`M` — never a silent wrong answer
       (Cardinal Rule 1).
   * - ``GreenOperator(-S, ...)``
     - *(refused)*
     - :math:`-S` is not invertible; the constructor raises
       ``NotInvertible`` at **construction**.

``is_invertible`` on a combination is now spelling-independent: a sum of
two or more terms is never invertible by itself, whatever operand order
spelled it, and the keystone faithfulness contract
(:ref:`capability-set-semantics`) still pins that ``.inverse()`` returns
exactly when ``is_invertible`` is ``True``.

.. dropdown:: First got wrong — the spelling chose the algorithm (until 2026-10-20)
   :color: muted

   From step 4 until package S2, a sum's inverse was keyed on its operand
   ORDER (the #261 canonical-ordering rule extended to the sum inverse):
   ``OperatorSum.is_invertible`` returned ``self.a.is_invertible``, the
   first-spelled term, and ``OperatorSum.inverse()`` built a
   ``GreenOperator`` whose preconditioner was the left-spine head, its
   other terms the gains (an ``A - S`` spelling arrived as
   ``ScaledOperator(-1, S)`` and was un-wrapped). The four edges were:
   ``L + C`` fused by the type test in ``StreamingOperator.__add__`` and
   inverted by the sweep; ``(L+C) - S`` built the physical Green; ``C + L``
   missed the fusion, built a Green with the collision head and diverged
   loudly; ``(-S) + A`` refused at construction, naming the canonical
   ordering. The predicate therefore read "leading-term-preconditionable at
   this operand order", not mathematical invertibility (verification spec
   §18.B), and that was accepted because #261 already made order
   load-bearing.

   It could not survive the formal linear combination: an unordered sum
   has no head. `[M]` the S2 census: no production path used it (the
   left-spine walk ran 52 times, all from tests; every production inverse
   of a loss goes through ``Splitting.implicit.inverse()``), and 5 tests
   pinned it, three on the message fragment "canonical ordering". Ruling
   50 moved the split into the constructor, ``GreenOperator(M, N)``, and
   ruling 52 removed the fusion that made ``L + C`` and ``C + L`` differ.


.. _green-strategy-mints:

The Strategy mints the Green
----------------------------

A sum of two or more terms has no inverse of its own
(:eq:`combination-inverse`), so the inverse of the within-group loss
exists only where a split has been named, and the split is the
Strategy's (ruling 65 of ``.claude/plans/numerics_untangling.md``). The
S\ :sub:`N` Strategy value is :class:`~orpheus.sn.splitting.Splitting`:
``Splitting.from_schedule`` labels each loss term implicit or explicit,
and the value derives :math:`M` (``implicit``) and the lagged pieces
:math:`N_i` (``explicit``) from the labels
(:ref:`sn-splitting-is-a-strategy-value`). The value then mints the
inverse it names, through one verb:
``Splitting.green(realisation, extra_gains=…)``
(:meth:`~orpheus.sn.splitting.Splitting.green`) returns
``GreenOperator(implicit, N, realisation=realisation)``, where :math:`N`
is the sum of the ``explicit`` pieces followed by ``extra_gains`` (the
multiplying-source production, which enters :math:`N` here and never
beside it). ``A.inverse()`` on the loss itself keeps refusing, and its
:class:`~orpheus.numerics.operator.NotInvertible` message names both
routes: ``GreenOperator(M, N)``, and in S\ :sub:`N` ``Splitting.green``.
An inverse with no Strategy in scope is therefore unspellable.

Every production solve of the within-group problem reads the Green the
Strategy minted, never a driver built beside it:

.. list-table:: Who reads the Strategy's Green in S\ :sub:`N` production
   :header-rows: 1
   :widths: 34 30 36

   * - Site
     - Green
     - Face read
   * - eigenvalue source-iteration inner
       (``SNSolver._solve_source_iteration``) and the fixed-source
       source-iteration solve (``_solve_fixed_source_si``)
     - ``_within_group_green_si``: a
       :class:`~orpheus.numerics.green_operator.Richardson` Green
     - ``green.driver(...)``
   * - eigenvalue Krylov inner (``SNSolver._solve_krylov``) and the
       fixed-source Krylov solve (``_solve_fixed_source_krylov``)
     - ``_within_group_green_krylov``: a
       :class:`~orpheus.numerics.green_operator.Krylov` Green
     - ``green.driver(...)``
   * - the eigenvalue finalize (:func:`~orpheus.sn.solver.solve_sn`)
     - the Green the last inner drove, read off the
       :class:`~orpheus.sn.solver.InnerSolve` record
     - ``green.step(q_F, ψ)`` (:ref:`sn-finalize-one-step`)
   * - the fixed-source finalize
     - the same Green as its solve
     - ``green.solution(q, ψ)``
   * - ``solve_sn_adjoint``, through
       :class:`~orpheus.numerics.iteration.KEigenvalue`
     - ``GreenOperator(M.H, N.H)``, the Green of the ADJOINT split
     - ``green.driver(...)``, built once by ``KEigenvalue``
   * - ``solve_sn_adjoint_fixed_source``
     - ``GreenOperator(M.H, N.H)``
     - ``green.driver(...)``

The adjoint entries mint the Green of the daggered split
:math:`(M^{*}, N^{*})` from the Strategy's own daggered parts; they do not
take ``.H`` of the forward Green (:ref:`green-adjoint-stays-base`).
Before package S4 the six driver sites each built ``SourceIteration(
seeded_inverse(M), N)`` or a ``KrylovAcceleration`` by hand, while the
``GreenOperator`` class had 1 production construction, its own ``dual``
(`[M]` the S4 census, 2026-10-20: 28 ``GreenOperator(`` call sites, 27 of
them in tests). The tree then held two realisations of one inverse, and
only the hand-wired one ran.

The S4 carve moved no production value: `[M]` the test-architect's
baseline of 26 production configurations (k, :math:`\phi`, :math:`\psi`
and the inner records) is bit-identical before and after
(``tests/gates/sn/solve/test_strategy_mints_green_sn.py``, docstring).
Because a value gate cannot see which route ran, the route is pinned by
counting spies (``test_g1a_every_production_driver_is_built_by_a_green``,
``test_g1b_every_finalize_is_one_step_of_the_green_the_last_inner_drove``)
and by an AST gate that the retired spellings are gone
(``test_g1c_the_retired_spellings_are_gone``).


.. _green-three-faces:

Three faces: the promise, the relaxation, the splitting map
-----------------------------------------------------------

The Green and its consumers keep two different contracts, and the Green
carries both as faces of one object (ruling 68), plus a third verb that
is not an inverse at all (ruling 69):

.. list-table:: The faces of ``GreenOperator``
   :header-rows: 1
   :widths: 24 40 36

   * - Face
     - Contract
     - Consumers
   * - :meth:`~orpheus.numerics.green_operator.GreenOperator.apply`,
       the PROMISE
     - the converged :math:`(M-N)^{-1}q`, with true relative residual below
       ``tol`` (:eq:`green-true-residual`), or a loud
       :class:`~orpheus.numerics.green_operator.ConvergenceFailure`
     - composites holding a Green, the invariant gates
   * - :meth:`~orpheus.numerics.green_operator.GreenOperator.driver`,
       the RELAXATION
     - a best-effort driver under a budget ``max_iter`` and a stopping
       tolerance ``tol``; its ``solve(q, initial_guess=…)`` returns the
       iterate and its :class:`~orpheus.numerics.convergence.IterationRecord`.
       A non-converged exit is legal and audible, never raised
     - the four within-group solve sites, the ``KEigenvalue`` inner, the
       adjoint fixed-source solve
   * - :meth:`~orpheus.numerics.green_operator.GreenOperator.step`,
       the SPLITTING MAP
     - :math:`G(\psi) = M^{-1}(q + N\psi)` applied ONCE, seeded by
       :math:`\psi`, through :math:`M`'s canonical seeded inverse
     - the two finalizes (reconstruction from a converged iterate)

The relaxation is constructed in one place: ``driver`` delegates to the
realisation, and ``apply``'s refinement loop is written over
``self.driver(...)``. So the promise is the relaxation run until the true
residual is met, and no second construction of the driver exists to
drift.

**Why the k-eigenvalue inner reads the relaxation and not the promise.**
``K = GreenOperator(M, N) @ F`` is the normal FORM of the k-eigenvalue
problem, the eigen-operator whose dominant eigenvalue is :math:`k`. Power
iteration does not need :math:`K` exactly at each outer step: the inner
solve need only be converged enough that the outer dominant-eigenvector
direction is accurate, and warm-starting each inner from the previous
outer iterate amortises its cost (inexact power iteration). Partial
convergence there is by design. Reading ``apply`` would either over-solve
every inner (the promise drives to ``tol`` on the equation residual) or
import the promise's hard raise into a loop that legitimately runs
budget-bounded solves. So :class:`~orpheus.numerics.iteration.KEigenvalue`
holds the Green, ``KEigenvalue(posing, green, …)``, and builds its inner
once as ``green.driver(max_iter=max_inner, tol=inner_tol)``. The normal
form is the concept; the warm-started relaxation of the same Green is its
production realization. ``KEigenvalue`` types the Green as
:class:`~orpheus.numerics.iteration.SupportsRelaxation`, a protocol with
one member, so that ``iteration`` does not import ``green_operator``,
which sits above it in the concept order.

**The splitting map is a verb of the split, not an inverse.** At a fixed
point :math:`G(\psi^{\star}) = \psi^{\star}`, so one application of
:math:`G` to a converged iterate returns it, which is the whole
reconstruction the finalizes need. ``step`` uses :math:`M`'s canonical
seeded inverse (``implicit_inverse``, built once by the constructor) and
never the realisation's preconditioner: after a moment-windowed solve the
un-windowed sweep turns the moment iterate back into per-ordinate
:math:`\psi`, and after a power iteration the finalize hands it the
CONVERGED eigenvalue's fission source, so the returned flux solves the
equation the caller is told about. Its right-hand side is the one
spelling :func:`~orpheus.numerics.iteration.lagged_source` with the
driven gains, the same body the source iteration evaluates each pass.
The seed carries the curvilinear :math:`\psi_{1/2}` seed and the
reflective partner trace. ``solution(q, ψ)`` is the iterate read back in
the split's own space: natively the iterate itself, under the moment
window one ``step`` (:ref:`green-realisations`); the representation
decides, never a flag beside it. Until package S4 the map was the free
function ``fixed_point_step(inverse, gains, q, ψ)``, which took the
inverse and the gains as two separate arguments a caller could pair
wrongly.


.. _green-realisations:

The realisation is a value: Richardson and Krylov
-------------------------------------------------

How :math:`(M - N)^{-1}` is iterated is a value the Strategy chooses and
hands to the Green (ruling 68): a
:class:`~orpheus.numerics.green_operator.Richardson` or a
:class:`~orpheus.numerics.green_operator.Krylov`, both in
:mod:`orpheus.numerics.green_operator`. A Richardson Green and a Krylov
Green of one split are the SAME operator (:ref:`green-operator`); the
realisation chooses the algorithm, and the two build the public drivers
:class:`~orpheus.numerics.iteration.SourceIteration` and
:class:`~orpheus.numerics.iteration.KrylovAcceleration`, so the Green
re-implements no iteration mathematics (taxonomy §11.2, Pattern 5).
``GreenOperator(M, N)`` with no realisation is ``Richardson()``, the plain
source iteration preconditioned by :math:`M^{-1}`.

**Richardson's preconditioner is derived, because it fixes the answer.**
The source iteration with a preconditioner :math:`P` is
:math:`\psi_{n+1} = P(q + N\psi_n)`. If :math:`P` is invertible its fixed
point satisfies :math:`P^{-1}\psi^{\star} = q + N\psi^{\star}`, that is

.. math::
   :label: green-richardson-fixed-point

   \psi^{\star} \;=\; \bigl(P^{-1} - N\bigr)^{-1} q ,
   \qquad\text{so}\qquad
   \psi^{\star} = (M - N)^{-1}q \ \ \forall q
   \iff P = M^{-1} .

.. (vv-status rationale) A two-line algebraic identity on the fixed point
   of the preconditioned source iteration; the S4 gates
   test_g9a_richardson_accepts_no_preconditioner and
   test_g9b_the_native_richardson_iterates_with_m_s_own_inverse pin the
   design consequence (no free P), not the identity.
.. vv-status: green-richardson-fixed-point documented

A Richardson preconditioner is therefore not a free parameter of the
iteration: it is the implicit part of a split, and a Richardson
realisation handed another split's inverse silently computes that
other split's Green. The design of S4 first took :math:`P` as a free
argument, ``Richardson(preconditioner, representation)``, and the
elegance review measured the consequence: `[M]` 2026-10-21,
``GreenOperator(D, N, realisation=Richardson(seeded_inverse(D2)))`` gave a
driver record reading "converged" after 30 iterations while the iterate
was 15.5 % off :math:`(D - N)^{-1}q`; ``apply`` raised
:class:`~orpheus.numerics.green_operator.ConvergenceFailure`, but the six
production sites read the driver face, which stays silent
(``scratch/numerics_untangling/packages/S4_review/elegance.md``, attack 1).
Ruling 70 removed the argument: the Green builds :math:`M^{-1}` once
(``implicit_inverse``), and ``Richardson(*, representation, corrector)``
derives :math:`P` from it through the representation, so a source
iteration over another split's inverse cannot be spelled. Every one of
the 3 production Richardson preconditioners before the ruling was already
a function of :math:`M` alone (the seeded inverse, the sweep, or the
windowed sweep), so nothing in production lost a capability.

**Krylov's preconditioner is free, because it does not.** GMRES with a
left preconditioner :math:`P` minimises :math:`\lVert P(q - (M-N)x)\rVert`
over a Krylov space of :math:`P(M - N)`. At convergence the preconditioned
residual vanishes, and for a nonsingular :math:`P`

.. math::
   :label: green-krylov-fixed-point

   P\bigl(q - (M - N)x^{\star}\bigr) = 0
   \;\Longrightarrow\;
   x^{\star} = (M - N)^{-1} q ,

.. (vv-status rationale) The left-preconditioned GMRES solution is
   independent of the nonsingular preconditioner; the gate that P moves
   the trajectory and never the converged fixed point is D4 of the
   Krylov-preconditioner suite (tests/gates/sn/solve/test_krylov_sweep_preconditioner.py).
.. vv-status: green-krylov-fixed-point documented

whatever :math:`P` is: the preconditioner moves the trajectory (the
Krylov space, hence the iteration count) and never the answer. So
``Krylov(P, restart=…)`` takes :math:`P` as the Strategy's free choice:
the sweep :math:`M^{-1}`, or the DSA-corrected
:math:`(I + \mathcal C)M^{-1}` (:ref:`sn-krylov-sweep-preconditioner`).
GMRES iterates in the split's own space, so a Krylov realisation's
representation is ``NATIVE``. ``restart`` is sized by the Strategy to
the full state the driver ravels (ERR-053).

**The representation is the space the iterate lives in.**
:class:`~orpheus.numerics.iteration.Representation` is a protocol of four
maps, all consequences of one choice, so no two of them can be paired
wrongly:

.. list-table:: The four maps of a ``Representation``
   :header-rows: 1
   :widths: 18 32 50

   * - Map
     - ``NATIVE`` (the split's own space)
     - :class:`~orpheus.sn.operators.windowing.AngularMomentWindow`
       (2-D Cartesian S\ :sub:`N`)
   * - ``precondition(M⁻¹)``
     - :math:`P = M^{-1}`
     - :math:`P = \mathcal P\,M^{-1}`, the fused
       :class:`~orpheus.sn.operators.windowing.WindowedSweep`
   * - ``rebind(g)``
     - the gain as it is
     - an angular lift re-bound on the moment iterate
       (``on_moment_domain()``); any other gain (the boundary gain, or
       its Gauss-Seidel upper part, which reads the trace) as it is
   * - ``analyse(state)``
     - the state as it is
     - :math:`\mathcal P` applied to a per-ordinate state (a seed)
   * - ``reconstruct(step, q, ψ̂)``
     - the iterate as it is
     - one splitting step, ``step(q, ψ̂)``, through :math:`M`'s own
       un-windowed inverse

The window holds exactly one datum, the analysis factor :math:`\mathcal P`
(a :class:`~orpheus.sn.operators.windowing.BulkAnalysisOperator` over the
scattering operator's own flux-analysis face), and derives the rest. Why
the windowed iteration computes the same Green: write
:math:`\tilde N` for the re-bound gains. Each angular lift reads the
per-ordinate flux only through the moments :math:`\mathcal P` keeps, and
the boundary gain reads the trace, which the windowed iterate carries
un-reduced, so :math:`\tilde N\mathcal P = N` on the states the
iteration produces. Then the windowed iterate is the analysis of the
native one, step for step: if :math:`\hat\psi_n = \mathcal P\psi_n`,

.. math::

   \hat\psi_{n+1}
   = \mathcal P M^{-1}\bigl(q + \tilde N\hat\psi_n\bigr)
   = \mathcal P M^{-1}\bigl(q + N\psi_n\bigr)
   = \mathcal P\psi_{n+1},

and at the fixed point :math:`M^{-1}(q + \tilde N\hat\psi^{\star}) =
M^{-1}(q + N\psi^{\star}) = \psi^{\star}`, which is ``reconstruct``: one
``step``. The solver chooses the representation at one site,
``_angular_representation`` in :mod:`orpheus.sn.solver`, on the genuine
condition ``is_cartesian and ndim == 2``; curvilinear meshes stay native
because the Morel–Montry Carlson seed reads the per-ordinate iterate at
:math:`\mu = -1`, which the moment tensor does not carry
(:ref:`sn-angular-windowing-geometry-restriction`). Before ruling 70 the
windowed sweep and the representation were chosen apart (the
``_maybe_window`` factory returned the sweep, and a bool beside it picked
the representation), so a windowed sweep with ``NATIVE``, or a plain
sweep with the window, could both be spelled (the elegance review,
attack 2).

**The driven gains are derived, scaled terms included.** A realisation
carries no gains. The Green reads :math:`N`'s terms
:math:`c_i w_i` in insertion order and drives
:math:`c_i\cdot\mathrm{rebind}(w_i)` (the word itself when
:math:`c_i = 1`), so every driven gain is a representative of a term of
:math:`N` by construction (the orchestrator's ruling on the S4 spec's Q3,
2026-10-20). The terms stay separate gains: ``GreenOperator(M, S + F)``
drives ``(S, F)`` and associates :math:`(q + S\psi) + F\psi`, exactly as
a hand-built ``SourceIteration(M⁻¹, S, F)`` does, where one combined
gain would associate :math:`q + (S\psi + F\psi)` (the S2 spec's finding
10). The rebind is applied to the WORD and the coefficient
outside it: a scaled lift :math:`c\,S` is an
:class:`~orpheus.numerics.operator.OperatorSum`, not an
:class:`~orpheus.transport.operators.angular_lift.AngularLift`, so
re-binding the scaled term misses the lift and the windowed driver feeds
a moment iterate to a per-ordinate binding. `[M]` the elegance review
(attack 5) found ``rebind(c * w)`` in the first S4 draft; `[M]`
2026-10-21, the archivist's probe on the 2-D heterogeneous
:math:`P_1` two-group fixture (``_cartesian_2d_p1_aniso_het_si``,
Jacobi): with :math:`N = 0.9\,S + N_{2n} + B` the windowed driver now
runs its 3 budgeted passes, where the review's probe raised
``TypeError`` from the scattering binding. Gated by
``test_g2d_a_scaled_angular_lift_is_rebound_then_scaled``.

.. important::

   **Gotcha — a Green's budget counts different units per
   realisation.** ``apply`` spends its ``max_iter``, the TOTAL budget, in driver passes,
   counting the cold start :math:`P q` as one. A Richardson pass is one
   splitting step; a Krylov pass is one GMRES restart cycle of up to
   ``restart`` Arnoldi steps. `[M]` 2026-10-21 (qa's
   ``S4_review/q1_paths.py`` re-run on the S4 worktree; 12-dof dense
   split, ``tol = 1e-12``): ``max_iter = 2`` with ``restart = 12``
   converged after 11 Arnoldi steps in one cycle, while the Richardson
   control at ``max_iter = 3`` ran 2 steps. A budget sized for one
   realisation is not a budget for the other. Giving the budget one unit
   is `#629 <https://github.com/deOliveira-R/ORPHEUS/issues/629>`_.

.. important::

   **Gotcha — the dual of a Strategy's Green is realised by the default
   Richardson.** ``G.dual()`` builds ``GreenOperator(M.dual(), N.dual())``
   with the same budget and no realisation: a Strategy's realisation is a
   choice for the forward split and is not transposed. So a Krylov Green
   whose split does not contract has a dual that diverges. `[M]`
   2026-10-21 (qa's ``S4_review/q1b_dual.py`` re-run on the S4 worktree:
   a 4×4 diagonal :math:`M` with a reversal-permutation :math:`N` scaled by
   6, :math:`\rho(M^{-1}N) = 1.90`): the forward Krylov Green converges to
   1.4e-15 relative, and ``dual().apply`` and ``apply_transpose`` both raise
   :class:`~orpheus.numerics.green_operator.ConvergenceFailure` (loudly,
   never a wrong answer). Gated by
   ``test_g9e_the_dual_of_a_strategy_green_carries_no_forward_realisation``.
   Transposing a Krylov realisation where its preconditioner transposes is
   `#629 <https://github.com/deOliveira-R/ORPHEUS/issues/629>`_.


.. _green-seed-law:

One seed law for a scaled inverse
---------------------------------

A one-term combination is invertible, :math:`(\alpha A)^{-1} =
\alpha^{-1}A^{-1}` (:eq:`combination-inverse`), and the inverse is a
first-class forward, the combination :math:`c\,w` with
:math:`c = \alpha^{-1}` and :math:`w = A^{-1}`. Its seed has one law
(ruling 67): ``initial_guess`` :math:`x_0` guesses the output
:math:`y = \alpha^{-1}A^{-1}b`, so the word's own output
:math:`A^{-1}b = \alpha y` is seeded with :math:`\alpha x_0`. The one-term
``apply`` forwards :math:`x_0 / c` to the word
(``OperatorSum._apply_scaled``), which is that law, and the transpose
face :math:`((\alpha A)^{-1})^{\mathsf T} = \alpha^{-1}(A^{-1})^{\mathsf T}`
obeys it too.

Before package S4 a driver never saw this law:
:func:`~orpheus.numerics.iteration.seeded_inverse` keyed on
wrap-delegate membership alone, so a scaled inverse fell into the
seedless branch and ``_SeededExactApply`` dropped the seed, while the same
object's direct ``apply`` honoured it. One object had two seed paths
that disagreed. Ruling 67 keys ``seeded_inverse`` on ONE capability
predicate, :func:`~orpheus.numerics.iteration.accepts_seed`: a
wrap-delegate member accepts a seed, a one-term combination accepts one
exactly when its word does, and nothing else does. A scaled inverse over a
seeded word is returned as it is, and its law holds in the driver slot as
on its direct face. `[M]` qa's probe (``S4_review/q1c_seed.log``,
2026-10-21):

.. list-table:: ``accepts_seed`` on the algebra-closed inverses
   :header-rows: 1
   :widths: 30 30 40

   * - Inverse
     - Type
     - ``accepts_seed``
   * - :math:`D^{-1}` (diagonal leaf)
     - ``InverseOperator``
     - ``True``
   * - :math:`(2D)^{-1}`, :math:`(-D)^{-1}`, :math:`(2\,D D)^{-1}`
     - ``OperatorSum`` (one term)
     - ``True``
   * - :math:`I^{-1}`
     - ``IdentityOperator``
     - ``False``
   * - :math:`(2I)^{-1}`
     - ``OperatorSum`` (one term)
     - ``False``, so ``seeded_inverse`` wraps it in ``_SeededExactApply``

So the sentence "an exact inverse has nothing to seed" is true of a
SEEDLESS word only. It is the whole story for the identity, a
permutation and a scaled identity; it is false for a scaled inverse over
the curvilinear sweep, whose Carlson closure reads the seed's
:math:`\mu = -1` flux, so a dropped seed there would be a wrong-fixed-point
defect rather than a rate change (`[R]`: no production site builds a
scaled inverse over the sweep today, which is why the old drop was
silent). The forward face, the transpose face
and the driver slot are each pinned by a gate with a reddening mutation
(``tests/gates/numerics/test_strategy_mints_green.py``,
``test_g4a`` to ``test_g4f``). `[M]` the S4 gate spec's finding F3: the
26-configuration production baseline is value-blind to the driver-slot
seed, which is why those gates record the seeds an exact inverse is handed
rather than comparing values.


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
physical :term:`scattering ratios <scattering ratio>` — because increment-stop delivers only
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


.. _green-adjoint-stays-base:

Why ``G.H`` stays the base
--------------------------

The Euclidean transpose of a Green is the Green of the dual split,
:math:`((M - N)^{-1})^{\mathsf T} = (M^{\mathsf T} - N^{\mathsf T})^{-1}`,
which :meth:`~orpheus.numerics.green_operator.GreenOperator.dual` builds
exactly. The Hilbert adjoint is the base's ``♯ ∘ dual ∘ ♭`` over that
dual, and S4 does not override it (ruling 66, until package S6), for two
reasons. First, the tempting override ``GreenOperator(self.inner.H)``
cannot be built: a Green needs a split, and the split cannot be read back
off ``(M - N).H``, which is a formal sum with no implicit part. Second,
the distributed form ``GreenOperator(M.H, N.H)`` is the honesty gate's
independent witness (``_sn_green`` in
``tests/gates/numerics/test_adjointability_is_honest.py``): an override of
that form was removed on 2026-10-18, because the witness then agreed with
``.H`` to exactly 0.0, a tautology. The law
:math:`((M - N)^{-1})^{*} = (M^{*} - N^{*})^{-1}` holds and is checked; it
is not a second route to ``.H``. That is also why the adjoint entries of
:ref:`green-strategy-mints` mint ``GreenOperator(M.H, N.H)`` from the
Strategy's daggered parts: the Green of the adjoint split is the
Strategy's value, not the forward Green's ``.H``. Ruling 53's distributed
adjoint is revisited at S6, where it holds everywhere. Gated by
``test_g6a_green_does_not_override_its_adjoint`` and
``test_g6b_no_production_green_is_built_on_daggered_parts_outside_the_strategy``.


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
The **transpose** axis is part of the back-half, as one law,
:math:`(A^{-1})^{\mathsf T} = (A^{\mathsf T})^{-1}` (transpose
:math:`A A^{-1} = I`): the mixin's ``apply_transpose`` is the forward's
transposed solve ``inner.solve_transpose``, and its ``is_adjointable`` is
:func:`~orpheus.numerics.operator.transpose_invertible` of the forward
(:attr:`~orpheus.numerics.operator.LinearOperator.is_transpose_invertible`,
a capability of the forward beside ``is_invertible`` and
``is_adjointable``). Read the other way, the inverse's own
``solve_transpose`` is the forward's ``apply_transpose`` and its
``is_transpose_invertible`` is the forward's ``is_adjointable``. The
mixin also passes the forward's ``is_metric_free_adjoint`` through, since
an operator that commutes with every metric (a constant multiplier
:math:`cI`) has an inverse that does too. ``.H`` itself stays the base's ``♯ ∘ dual ∘ ♭``, built when the
transpose exists and the adjoint is posed
(:ref:`operator-adjoint-two-questions`). Two siblings override the law:
:class:`~orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator`
(the ``trans=1`` backsolve on its own LU factors) and ``GreenOperator``
(its dual is the Green of the dual split, since a sum has no transposed
solve of its own; its ``.H`` is the base's).


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
       hand-built ``SourceIteration(sweep, S)`` (``array_equal``; the
       lagged part's terms ride as separate gains), and the edges of the
       split table above end-to-end: the plain sum refuses an inverse in
       both spelling orders, and the collision-preconditioned split
       ``GreenOperator(C, -L)`` raises ``ConvergenceFailure``.
   * - ``test_strategy_mints_green.py`` (numerics) and
       ``test_strategy_mints_green_sn.py`` (S\ :sub:`N`), package S4
     - The ROUTE: counting spies that every production driver is built by
       ``GreenOperator.driver`` and every finalize is one ``step`` of the
       Green the last inner drove, and an AST gate that ``fixed_point_step``
       and the hand-built drivers are gone (G1). The driven gains are
       :math:`N`'s terms through the representation, a scaled lift
       re-bound then scaled (G2). ``step`` is one application of the
       splitting map, seeded by the iterate, and reproduces a converged one
       (G3). The seed law of a scaled inverse on the forward face, the
       transpose face and the driver slot, keyed by ``accepts_seed`` (G4).
       The promise over a Krylov realisation, converged or raised, and
       written over the relaxation (G5). ``G.H`` not overridden (G6).
       ``KEigenvalue``'s inner built through the Green's driver (G7).
       Richardson accepts no preconditioner and derives it from :math:`M`;
       a Krylov Green's preconditioner is the Strategy's (G9).

The config discipline is Mode-9: the L1 slab is **heterogeneous, ≥2-group,
vacuum** — a reflective isotropic box would null the
streaming↔scattering redistribution the Neumann series expands, and 1-group
is blind to a scattering-matrix transpose (Mode 6). The teeth are
**14 mutations verified** (2026-07-02, on the pre-S2 left-spine Green;
the test module records that they re-run at the S2 review): 12 bite their
named gates
(sign/swap/flatten/tol/increment/seed/order/…), and 2 are designed-green
controls — the always-wrap ``_negated`` no-op (a helper retired at S2;
it proved the gain-unwrap was pure deforestation) and the **M-GRN-SEED blindness proof** (a dropped
driver-start reddens Green's own Mode-11 seed spy while the landed step-3
spy and every value gate stay green, proving the step-3 spy is structurally
blind to Green's driver-start threading). The pyright ratchet held exactly
at baseline; the Sphinx build stayed ``-W`` clean.


.. _matrix-inverse-operator:

The materialising functor and the dense direct inverse
======================================================

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

**The return is a dense** :class:`numpy.ndarray`, and ``as_matrix`` is
realised by the apply-to-basis loop above. Three classes override it, each
reading a matrix it already holds rather than computing one from structure:
:class:`~orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator`
runs a single batched LU back-solve (below),
:class:`~orpheus.numerics.assembled_operator.SparseAssembledOperator`
densifies its sparse matrix, and
:class:`~orpheus.numerics.coupled_system.CoupledOperator` probes over typed
block-field basis elements. The *sparse* realization is not an override of
``as_matrix``: it is a second algorithm for the same matrix, structural
assembly, which lives downstream of the algebra as
:func:`~orpheus.numerics.assembled_operator.assemble` (sparse) and
:func:`~orpheus.numerics.assembled_operator.dense` (densified), and the
consumer that pays for a densification chooses between the two
(:ref:`assembly-matrix-and-algorithms`). No octant-batched 3-D streaming
assembly exists; the SN assembly mode emits Cartesian per-ordinate blocks
(:ref:`loss-rep-three-modes`).


``MatrixInverseOperator`` — the dense direct inverse
----------------------------------------------------

``A.inverse()`` is a **structure-keyed factory** (taxonomy §13): the
concrete subclass names the mathematical *object* the inverse is, never the
algorithm. A schedule-triangular forward returns the direct-substitution
:class:`~orpheus.sn.operators.sweep_operator.SweepOperator`; a sum of two
or more terms has no inverse of its own, and the preconditioned-splitting
:class:`~orpheus.numerics.green_operator.GreenOperator` is built from its
named split, ``GreenOperator(M, N)``; a value-bearing
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
whose only spelled split has a **non-invertible implicit part**: the
motivating form is the SN :math:`(-S)+(L+C)`, which as a sum has no
inverse of its own (:eq:`combination-inverse`), and whose split with
:math:`M = -S` :class:`~orpheus.numerics.green_operator.GreenOperator`
**refuses at construction**
(:class:`~orpheus.numerics.operator.NotInvertible`), yet which
materializes to a perfectly invertible matrix. Because the real :math:`(-S)+(L+C)` is a
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
  rather than through an iterative inverse: the loss
  :math:`C - K_{\rm iso}` is a combination with no inverse of its own, and
  its iterative inverse would be the Green of a named split,
  ``GreenOperator(C, K_iso)`` — the direct-realization strategy encoded as
  a type. With the operator
  spelling live, ``direct_eigenvalue`` has **zero production consumers**,
  retained as the ``(A, F)``-posed sibling engine and the RQI test oracle.
  See :ref:`direct-eigensolve-solve`.
* **CP :math:`[P]` — the next production method in line.** The
  collision-probability matrix is dense **by construction**
  (:doc:`/theory/methods/collision_probability`, §14b), so a CP ``.inverse()`` realized as a
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
:meth:`~orpheus.transport.operators.isotropic_transfer.IsotropicScattering.dense_per_material`
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
   * - ``tests/gates/numerics/test_matrix_inverse_operator.py``
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

The assembly functor — structural sparse emission
=================================================

An operator's matrix is fixed by its spaces: the domain's and the
codomain's coordinates, the axes in their order, decide which row and
which column every entry occupies, and
:meth:`~orpheus.numerics.operator.LinearOperator.as_matrix` computes it by
probing (:eq:`matrix-functor-out`). An operator that knows its own
:math:`(\text{row},\text{col},\text{value})` footprint (a finite-difference
stencil, a cell-block-diagonal transfer, a multiplier's diagonal) can
compute the same matrix faster, by scattering its coefficients into a
sparse matrix instead of applying itself :math:`n` times. That second
algorithm is **assembly**, and this section is its operator-algebra view:
where the algorithm lives, how it threads through the composites, and how
the two algorithms are held to one matrix. The per-method cell mathematics
(the SN closure walk, the diffusion finite-difference stencil) is developed
in :doc:`/theory/methods/sn/loss_representation`
(:ref:`loss-rep-three-modes`) and :doc:`/theory/methods/diffusion_1d`.

Everything named here lives in :mod:`orpheus.numerics.assembled_operator`
except one property: the structural flag
:attr:`~orpheus.numerics.operator.LinearOperator.is_assemblable`, which
:mod:`orpheus.numerics.operator` states and which is the only thing the
operator algebra knows about assembly. The design is package S3 of the
numerics untangling (rulings 24 and 60–62 of
``.claude/plans/numerics_untangling.md``, 2026-10-20).

.. _assembly-matrix-and-algorithms:

The matrix and its two algorithms
---------------------------------

The user's framing at the ruling (ruling 60) was: *"if space knows its
axis, and LinearOperator knows its space, LinearOperators should have no
problem assembling the associated matrix"*. It resolves into two separate
facts. **The spaces fix where each entry sits**: the flat layout of the
domain and the codomain is a property of the spaces, so the matrix is a
property of the operator and its spaces, defined for every operator that
has them. **The algorithm fixes how each entry is computed**: probing or
assembly, two procedures for the one matrix. The code puts each fact in
its own place:

.. list-table::
   :header-rows: 1
   :widths: 22 26 26 26

   * - Algorithm
     - Spelling
     - Domain
     - Cost
   * - **probing** (apply to the unit basis)
     - :meth:`~orpheus.numerics.operator.LinearOperator.as_matrix`
     - **total**: every operator with spaces has a matrix
     - dense; :math:`n` operator applications
   * - **assembly** (structural emission)
     - :func:`~orpheus.numerics.assembled_operator.assemble` (sparse) and
       :func:`~orpheus.numerics.assembled_operator.dense` (densified)
     - **partial**: only where every part emits
     - sparse; :math:`O(\mathrm{nnz})` scatter

``as_matrix`` is the representation and is realised by probing alone: it
never consults ``is_assemblable``. Three classes override it, and each
override reads a matrix the object already HOLDS rather than computing a
new one by structure:
:class:`~orpheus.numerics.assembled_operator.SparseAssembledOperator`
densifies its CSR matrix;
:class:`~orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator`
runs one batched LU back-solve against the identity;
:class:`~orpheus.numerics.coupled_system.CoupledOperator` probes over typed
:class:`~orpheus.numerics.coupled_system.CoupledField` basis elements,
because a typed block refuses the bare flat columns the base probe feeds
it. An operator that only KNOWS its stencil does not override
``as_matrix``. Its structural densification is ``dense(op)``, which
assembles when ``op.is_assemblable`` holds and otherwise falls back to
``op.as_matrix(...)``, under the same basis-shape resolution and the same
:class:`~orpheus.numerics.operator.MatrixTooLarge` size gate (one helper,
``_check_dimension`` in :mod:`orpheus.numerics.operator`, gates both). The
consumer that pays for a densification chooses the algorithm; in
production that consumer is ``MatrixInverseOperator``
(:ref:`assembly-why-matrix-inverse-chooses-dense`).

The separation has a verification consequence. Because ``as_matrix`` never
routes through assembly, ``op.as_matrix()`` and ``dense(op)`` on an
assemblable operator are two independent algorithms, and comparing them
is a comparison (:ref:`assembly-verification`), not assembly compared with
itself.

The carrier
-----------

Assembly lands in
:class:`~orpheus.numerics.assembled_operator.SparseAssembledOperator`, a
thin :mod:`scipy.sparse` wrapper conforming to
:class:`~orpheus.numerics.operator.LinearOperator`: ``apply`` is the CSR
matvec on a flat vector, ``apply_transpose`` the exact CSR transpose
(``is_adjointable`` is ``True``), and ``is_invertible`` is the base
``False`` (an inverse is an explicit consumer choice, such as a
``MatrixInverseOperator`` over it). It is **not** a COO-builder type with
an algebra of its own: that would restate every law of the operator
algebra on triplet buffers, one layer down. Its constructor converts any
scipy sparse input to CSR once, and scipy's ``COO → CSR`` conversion SUMS
repeated :math:`(\text{row},\text{col})` entries, the finite-element
scatter semantics, so an emitter may scatter per-cell and per-face
contributions freely. Its ``domain`` and ``codomain`` are passthrough
metadata from the operator it assembles. The functor is closed on the
carrier: ``assemble(S) is S`` for an assembled operator ``S``.

The assembly surface
--------------------

The assembly axis has the same three layers as the inverse and adjoint
axes (:ref:`capability-set-semantics`), plus a refusal, but only its first
layer lives in the algebra:

#. a **predicate**,
   :attr:`~orpheus.numerics.operator.LinearOperator.is_assemblable`
   (base default ``False``; an operator emits only by explicit override).
   Every composite derives it from its parts: a combination emits iff every
   term does, a word iff every factor does, a fused atom iff its expansion
   does, and each domain composite iff its parts do. Examples of operators
   that emit nothing: :class:`~orpheus.numerics.operator.IdentityOperator`,
   :class:`~orpheus.numerics.operator.ZeroOperator`, the numerics
   :class:`~orpheus.numerics.elementary_operators.DiagonalOperator`,
   :class:`~orpheus.numerics.operator.TensorProductOperator`, and the SN
   sweep atoms, whose expansions do not emit;
#. an **emission verb**, ``assembled_matrix(emit)``, defined on each
   emitter and returning a RAW scipy sparse matrix in the flat layout of
   the operator's spaces. ``emit`` is the functor's own raw emission,
   handed to the domain composites that recurse into their parts; a true
   leaf ignores it. The narrowing Protocol that declares the verb is
   :class:`~orpheus.numerics.assembled_operator.SupportsAssembly`, and the
   checked bridge that narrows to it is the ``TypeGuard``
   :func:`~orpheus.numerics.assembled_operator.assemblable` (never
   ``isinstance`` the Protocol: narrow with the bridge);
#. an **operator-returning map**, the functor
   :func:`~orpheus.numerics.assembled_operator.assemble`, which returns
   :math:`[A]` wrapped once in a ``SparseAssembledOperator`` carrying
   ``A``'s ends.

The refusal is
:class:`~orpheus.numerics.assembled_operator.MissingAssembly`, a
``TypeError`` and the assembly-axis sibling of
:class:`~orpheus.numerics.operator.NotInvertible` and
:class:`~orpheus.numerics.operator.MissingAdjoint`. The functor raises it
when it reaches a part that emits nothing, and its message names that
part's type and the law that needed it: ``assemble(a + IdentityOperator())``
refuses with ``[Σ cᵢwᵢ] = Σ cᵢ[wᵢ]: IdentityOperator emits no assembly``.
A domain composite whose own flag is ``False`` refuses at itself, since
its flag already says that some part cannot emit. ``dense`` does not
refuse: an operator whose flag is ``False`` is densified by probing.

The emitters are the diffusion
:class:`~orpheus.diffusion.operators.LeakageOperator` and
:class:`~orpheus.diffusion.operators.DiffusionBoundaryOperator`, the
transport
:class:`~orpheus.transport.operators.isotropic_transfer.IsotropicTransfer`
(inherited by the isotropic scattering and :math:`(n,2n)` operators) and
:class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`,
the carrier itself, and the three domain composites of the next
subsection. The SN per-ordinate blocks are assembled outside the functor,
by :func:`~orpheus.sn.loss_representation.assembly.assemble_ordinate_blocks`,
which builds its carriers directly
(:ref:`loss-rep-three-modes`).

The homomorphism laws — one functor downstream of the algebra
-------------------------------------------------------------

Assembly is a homomorphism out of the operator algebra into the algebra of
sparse matrices. The free algebra of :ref:`composition-algebra` has three
composites, a formal linear combination of terms
(:class:`~orpheus.numerics.operator.OperatorSum`), a flattened word of
factors (:class:`~orpheus.numerics.operator.OperatorProduct`), and a fused
atom (:class:`~orpheus.numerics.operator.FusedAtom`, an operator its owner
built as one term, carrying an ``expansion()`` into the other two), and
the functor carries one law per composite:

.. math::
   :label: matrix-functor-homomorphism

   \Bigl[\,\sum_i c_i\,w_i\Bigr] = \sum_i c_i\,[w_i], \qquad
   [A_1 A_2 \cdots A_n] = [A_1]\,[A_2]\cdots[A_n], \qquad
   [\,a\,] = [\operatorname{expansion}(a)],

where :math:`[\cdot]` is the sparse matrix in the flat layout of the
operator's spaces, :math:`c_i` are the real coefficients of a combination,
:math:`w_i` its terms, :math:`A_1 \cdots A_n` the factors of a word and
:math:`a` a fused atom. The familiar binary laws are special cases:
:math:`[A+B] = [A] + [B]` is a two-term combination with unit
coefficients, :math:`[\alpha L] = \alpha[L]` the one-term combination (a
scaled operator IS a one-term combination since package S2), and
:math:`[AB] = [A][B]` a two-letter word.

The laws are realised once, in the private ``_emit`` of
:mod:`orpheus.numerics.assembled_operator`, as a ``match`` on the three
composites followed by the leaf arm:

* the **combination arm** folds the terms' emissions LEFT, in the order
  the terms were first inserted, with the same helper ``_left_fold`` that
  folds the combination's ``apply`` (out of place, because a sparse
  matrix's in-place operators are not its out-of-place ones). The matrix
  of a sum is therefore rounded in the same order as its action;
* the **word arm** multiplies the factors' emissions left to right with
  CSR ``@``;
* the **fused-atom arm** assembles the atom's expansion, so a capability
  the expansion gains reaches the atom with no edit to the atom;
* the **leaf arm** refuses a part whose flag is ``False``, and otherwise
  calls the part's ``assembled_matrix(emit)`` and converts the result to
  CSR before any law touches it
  (:ref:`assembly-why-csr-at-the-leaf`).

``assemble(op)`` returns ``op`` itself when ``op`` is already a
``SparseAssembledOperator``, and otherwise wraps ``_emit(op)`` ONCE. No
emitter constructs the carrier: a leaf returns its raw matrix, and only
the functor wraps (ruling 62).

**The domain composites** (ruling 61) are three operators that are
composites in their own domain but atoms to the free algebra:
:class:`~orpheus.numerics.coupled_system.CoupledOperator` (a block grid
over a coupled space), :class:`~orpheus.transport.operators.lift.BulkLift`
(a bulk operator lifted onto bulk :math:`\oplus` trace) and
:class:`~orpheus.numerics.flat_operator.FlattenedOperator` (a typed
operator viewed on flat vectors). Each keeps its assembly law in its own
module and recurses through the ``emit`` it is handed:

* ``CoupledOperator`` scatters each present block's emission at its
  :math:`(\text{row}_i,\text{col}_j)` offset with
  :func:`scipy.sparse.block_array`, refusing with the offending block's
  index when one cannot emit (:ref:`coupled-block-operator`);
* ``BulkLift`` re-indexes its inner operator's emission into the
  composite flat layout (the bulk C-ravel followed by the trace buffer,
  zero trace rows and columns) with
  :func:`~orpheus.transport.operators.lift.embed_bulk_assembly`, which maps
  a raw matrix to a raw matrix and refuses a shape mismatch;
* ``FlattenedOperator`` passes its inner operator's emission through
  unchanged, because the typed operator's flat layout IS the flattened
  view's, after checking the emission's dimension against the template's
  flat size (a mismatch is an emitter drifting from the carrier's
  ``to_flat`` layout, and raises ``ValueError``). The carrier that
  ``assemble`` wraps around it carries the flattened view's ends, which
  are its inner operator's spaces.

No tensor-product law exists:
:class:`~orpheus.numerics.operator.TensorProductOperator` emits nothing,
so an expression containing one refuses with ``MissingAssembly``. The law
would be :math:`[A\otimes B] = [A]\otimes[B]`, a
:func:`scipy.sparse.kron`; no assembled operator contains a tensor-product
factor (the diffusion loss and the SN :math:`L + C` trees have none).

Why assembly is a homomorphism out of the algebra
-------------------------------------------------

**The free algebra must not know matrices.** The operator algebra is a
free algebra over its atoms: a combination and a word are formal objects,
and their meaning is given by the homomorphisms out of them. ``apply``
(the action on a vector) is one, the dual is another, and assembly into
sparse matrices is a third. If each composite carried its own law as a
method, the algebra would import its target, the sparse carrier, and
every new representation would add one method to every composite: the
laws would be written once per composite instead of once per
representation. A functor downstream of the algebra states each law once,
as an arm on the composite's structure, and leaves the algebra closed
under its own operations.

**The import cycle it removed.** While the laws lived in the composites'
bodies, :mod:`orpheus.numerics.operator` imported
``SparseAssembledOperator`` from :mod:`orpheus.numerics.assembled_operator`
(once under ``TYPE_CHECKING`` and lazily inside the sum's and the word's
assembly bodies), while ``assembled_operator`` imports ``LinearOperator``
from ``operator``: a cycle between the two modules. The concept-order gate,
``tests/gates/numerics/test_concept_order.py``, reads every ``import``
statement of the numerics package (lazy and ``TYPE_CHECKING`` imports
included) against the tier list ``TIERS``, where ``operator`` stands alone
in the lowest tier; the edge ``("assembled_operator", "operator")`` was
listed in its ``PENDING`` table as owed to S3. Today ``operator.py``
imports nothing from ``assembled_operator`` (its mentions of it are
docstring cross-references, which are not imports), the ``PENDING`` entry
is gone, and the gate reads the edge as a violation if it returns.

**The emitters do not import the functor.** A domain composite recurses
through the ``emit`` argument rather than calling ``assemble``. This is
what lets :class:`~orpheus.numerics.flat_operator.FlattenedOperator`
recurse at all: ``flat_operator`` and ``assembled_operator`` share one
tier of ``TIERS``, and a same-tier import is a violation of the gate. The
same argument keeps the emitters' modules free of the assembly module: of
them only ``coupled_system`` (one tier up) imports it, and only for
``MissingAssembly``.

.. _assembly-why-csr-at-the-leaf:

Why the CSR canonicalisation sits at the leaf arm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An emitter usually returns COO triplets with repeated
:math:`(\text{row},\text{col})` entries (a shared face contributes to the
same diagonal from both cells). The leaf arm converts each emission to CSR
before the combination arm multiplies it by a coefficient, so duplicates
are summed first and scaled second: :math:`c\,(a + b)`, not
:math:`c\,a + c\,b`. In floating point the two differ, and the first is
the rounding the carrier fixed when every leaf still returned a
``SparseAssembledOperator`` (whose constructor sums duplicates) and the sum
law read its CSR ``.matrix``. Converting at the leaf keeps every assembled
matrix bit-identical to what it was before the functor moved; converting
in ``assemble`` alone, after the laws, would scale unsummed duplicates. The
row ``test_g1g_a_leaf_emission_is_summed_before_its_coefficient`` in
``tests/gates/numerics/test_assembly_functor.py`` pins the order.

.. _assembly-why-matrix-inverse-chooses-dense:

Why ``MatrixInverseOperator`` is the consumer that chooses ``dense``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:class:`~orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator`
LU-factors a materialised matrix once, at construction, and it is the one
production consumer that densifies an operator which emits: the diffusion
resolvent ``MatrixInverseOperator(FlattenedOperator(A, template))``
(``orpheus/diffusion/solver.py``), the homogeneous solver's inverse of its
loss combination (``orpheus/homogeneous/solver.py``), and
``CoupledOperator.inverse`` (its dense direct route, on a square grid). It pays for the
densification, so it chooses the algorithm: ``dense(inner)``, which
scatters :math:`O(\mathrm{nnz})` coefficients where ``as_matrix`` would
apply the operator :math:`n` times. Its constructor is the one production
call of ``dense``; ``matrix_inverse_operator`` therefore imports
``assembled_operator`` and sits alone in its own tier of ``TIERS``, above
the assembly module.

The route also decides the matrix's bits. `[M]` 2026-10-20, the S3
verification spec's capture over ``tests/gates/{numerics,diffusion,homogeneous}``
(``-m "not slow"``), recording the SHA-256 of every matrix handed to
``lu_factor``: 176 matrices in 139 tests. A simulation of ruling 60 on
the pre-S3 tree, with ``MatrixInverseOperator`` assembling, changed 0 of
the 139 tests' digests; the same simulation with ``MatrixInverseOperator``
probing changed 9 of 139, because probing and assembly sum the stencil in
different orders. So ``dense`` is not only faster: it keeps the matrices
the diffusion and homogeneous solvers factor exactly as they were.

Alternatives that do not work
-----------------------------

Three other placements were put to the user at the S3 design checkpoint and
rejected there (2026-10-20, rulings 60 and 61). Each fails for a structural
reason, not a cost:

* **Moving** ``as_matrix`` **downstream, beside the functor.** The matrix
  is a property of the operator and its spaces, and probing needs only
  ``apply`` and the spaces, both of which the algebra already has. Moved
  downstream, the total representation would depend on the partial
  algorithm's module, every consumer of a dense matrix (the pencil, the
  diffusion boundary operator's per-law blocks, the derivations) would
  import the assembly layer although none of them assembles, and the
  representation would read as if it depended on an algorithm. The
  ruling keeps the representation with the spaces and the algorithm
  downstream.
* **A registration hook**: ``operator.py`` keeps an ``as_matrix`` that
  delegates to a hook, and ``assembled_operator`` registers the laws when it
  is imported. This hides the dependence instead of removing it. The
  algebra's behaviour would become runtime state set by import order: the
  same ``as_matrix`` call would probe or assemble depending on whether some
  other module had imported ``assembled_operator`` first. And the
  concept-order gate, which reads import statements, could not see the
  edge it exists to forbid.
* **All the domain composites' laws in** ``assembled_operator``. The block
  law, the lift embedding and the flat pass-through each need their
  composite's class. ``coupled_system`` sits one tier ABOVE
  ``assembled_operator`` in ``TIERS``, ``flat_operator`` shares its tier,
  and :class:`~orpheus.transport.operators.lift.BulkLift` lives in
  ``transport/``, layer L2 of the architecture map, above ``numerics/``
  (layer L1; :ref:`architecture-layering`). The numerics module would have
  to import upward, sideways and out of its layer. Ruling 61 is the
  structural answer: each domain composite is an atom to the free algebra
  and keeps its law in its own module, recursing through ``emit``.

.. _assembly-verification:

Verification — two algorithms, one matrix
-----------------------------------------

**The laws.** :eq:`matrix-functor-homomorphism` is verified at two levels.
``tests/gates/numerics/test_assembled_operator.py`` pins the binary laws
(sum, product, one-term scaling, difference) on exact dyadic leaves, where
every association gives the same bits. Those rows cannot see the fold
ORDER, so ``tests/gates/numerics/test_assembly_functor.py`` adds the arms on
a seeded non-dyadic fixture where the left fold of :math:`A - B + 0.3\,C`
differs from the other associations (`[M]` 2026-10-20, the S3 spec: 2
entries differ), the word order on integer-valued leaves, the fused arm,
nesting through every arm, the domain composites' recursion through
``emit``, the refusal's part and law, and the single wrap per functor call.

**Probing against assembly.** Because ``as_matrix`` is the probe and
``dense`` the assembly, a gate comparing them compares two algorithms. One
permanent pin per family:

* diffusion, ``test_g3_probed_equals_assembled_pin``
  (``tests/gates/diffusion/test_operators.py``): ``dense(flat)`` against
  ``flat.as_matrix()`` for the flattened diffusion loss;
* SN, ``test_g3_dd_slab_probed_column_pin``
  (``tests/gates/sn/sweep/test_assembly_mode.py``): the diamond-difference
  slab block from ``assemble_ordinate_blocks`` against columns probed by
  hand through ``A.apply``.

The two algorithms share an upstream by design: both read the same
coefficient sources (for diffusion the interior conductances, the boundary
closures, the multiplier's coefficient and the energy kernels), which is
the one-source discipline of :ref:`loss-rep-assembly-verification`. A
probed-against-assembled pin is therefore blind to an error in a
coefficient source. The absolute references cover that:
``test_assembled_loss_matches_hand_posed_stencil`` compares BOTH
``flat.as_matrix()`` and ``dense(flat)`` with the hand-posed stencil, and
the one-source teeth move all three consumption modes together.

The first production consumer — the diffusion resolvent
-------------------------------------------------------

The exact inner solve of the diffusion eigenvalue iteration is
``MatrixInverseOperator(FlattenedOperator(A, template))``: one eager LU at
construction, one back-substitution per outer iteration
(:doc:`/theory/methods/diffusion_1d`). Every diffusion loss leaf emits, so
``dense`` assembles the flattened loss and the LU factors the ASSEMBLED
matrix; the whole diffusion suite (its eigenvalue, trace and balance gates)
is the assembled path's regression net.

A value gate cannot tell which route ran, because on the production
fixtures the two routes agree bit for bit. The route is pinned by a
sentinel (vv-principles Mode 11, a green gate that never executes the
rewired line): ``test_resolvent_materializes_through_assembly`` spies on
``LinearOperator.as_matrix`` and on the functor, constructs the resolvent,
and asserts that the functor was called on the flattened loss and that no
``as_matrix`` probe fired at the composite flat dimension. Probes at the
law dimension :math:`n_g` are expected: the boundary emitter extracts each
realised boundary law's :math:`n_g\times n_g` block through the law's own
``as_matrix``, which is the one-source discipline (the law has no stencil
of its own), not a route escape.

`[M]` 2026-10-20, the S3 verification spec's route discriminator: on the
heterogeneous diffusion fixtures, probing and assembly agree to
:math:`\max|\Delta| = 0.0` in 3 of 3 configurations with the collision
built from :math:`\Sigma_t`; built from :math:`\Sigma_a` (a valid operator,
the fixture of the stencil-mutation tooth) they differ by 1 ULP
(relative :math:`7.5\times10^{-17}`) at 2 entries in 3 of 3
configurations. The bit-agreement is therefore a property of the fixture's
values, not of the two algorithms.

.. dropdown:: First got wrong — ``as_matrix`` delegating to assembly (2026-07-04 until 2026-10-20)

   When the assembly axis landed (stencil-assembly Phase 2b), each composite
   carried its law in an ``assemble()`` method, and ``as_matrix`` delegated
   to ``self.assemble().as_matrix(...)`` whenever the operator was
   assemblable. The representation then changed algorithm with the
   operator, so a probed-against-assembled gate would have compared assembly
   with itself; the probing loop had to be kept under a separate private
   name, ``_as_matrix_by_probing``, for the gates to force, and the Mode-11
   sentinel spied on that private name. Two concepts, the matrix and the
   algorithm that computes it, were welded into one method. Ruling 60
   separated them: ``as_matrix`` is the probe, ``dense`` the assembly, the
   consumer chooses, and the private hook retired because ``as_matrix``
   itself is the independent algorithm.


.. _operator-inverse-family-history:

Development history
===================

A reverse-chronological (latest first) changelog of the architectural
milestones of the inverse family. The table opens at package S4; the
family's earlier milestones (the #226 taxonomy steps, the wrap-delegate
mixin, the named split of package S2, the assembly functor of package
S3) are rows of :ref:`the operator algebra's history
<operator-algebra-history>`. Entries marked *(in development)* live on an
unmerged feature branch and have no landed merge-to-``main`` hash yet;
trust ``git`` over this table for merge status.

.. list-table::
   :header-rows: 1
   :widths: 10 50 12 28

   * - When
     - Architectural milestone
     - Issue
     - Where
   * - 2026-10-21
     - **The Strategy mints the Green, and the Green owns its
       realisation** (package S4 of the numerics untangling; rulings
       65–70). ``Splitting.green(realisation, extra_gains=…)`` mints
       ``GreenOperator(M, N, realisation=…)``, and every production
       within-group solve, the ``KEigenvalue`` inner and the adjoint
       entries read that Green instead of a hand-built
       ``SourceIteration(seeded_inverse(M), N)`` or ``KrylovAcceleration``.
       The Green has three faces: ``apply`` the promise, ``driver`` the
       relaxation (constructed in one place, and ``apply`` is written over
       it), ``step`` the splitting map, which retired the free function
       ``fixed_point_step``. The realisation is a value:
       ``Richardson(*, representation, corrector)`` derives its
       preconditioner from :math:`M^{-1}` through the ``Representation``
       (four maps: ``precondition``, ``rebind``, ``analyse``,
       ``reconstruct``; the 2-D window is ``AngularMomentWindow``, which
       replaced ``_maybe_window``), and ``Krylov(P, restart)`` keeps a free
       left preconditioner. The driven gains are derived from :math:`N`'s
       terms, a scaled lift re-bound then scaled. ``KEigenvalue(posing,
       implicit, explicit)`` became ``KEigenvalue(posing, green)``.
       ``seeded_inverse`` keys on one predicate, ``accepts_seed``, so a
       scaled inverse obeys one seed law in the driver slot and on its
       direct face. ``G.H`` stays the base ``♯ ∘ dual ∘ ♭``. Pinned by
       ``tests/gates/numerics/test_strategy_mints_green.py`` and
       ``tests/gates/sn/solve/test_strategy_mints_green_sn.py``. See
       :ref:`green-strategy-mints`.
     - —
     - ``[LANDED fb034da3]``
