---
name: operator-inverse-taxonomy-frames
description: TAXONOMY verdict for the #226 inverse-as-operator carve — substrate-vs-operator layer separation, the four inverse strategies, the dagger-partial-inverse category, and the layer-confusion smell (vestigial-forward resolvent). Operator-algebra spine.
metadata:
  type: project
---

# Operator-inverse taxonomy (#226 carve) — frame verdict

DESIGN session (no code): stress-test the user's substrate/operator-views
hypothesis for `A.inverse()`. Branch-direct reads (L-005), not Nexus-on-branch.
The hypothesis was ~85% already built + correct; two slogans broken.

## The load-bearing structural fact: TWO LAYERS, do not conflate

- **SUBSTRATE layer** = `LossRepresentation` (`sn/loss_representation/__init__.py:229`).
  ONE object, MANY action-views — `loss_action`(=`(L+C)ψ`, :253),
  `loss_action_transpose`(=`(L+C)ᵀφ`, :294), `sweep`(=`(L+C)⁻¹q`, :241/:500),
  `sweep_transpose`(=`(L+C)⁻ᵀr`, #280). It is NOT a morphism, so "three actions
  of one object" (the L21 claim) violates NOTHING here. This is where
  "one object, many views" is CORRECT.
- **OPERATOR layer** = `LinearOperator` (`numerics/operator.py:285`). Each is
  ONE morphism, ONE promise (`apply`). `.H`/`.inverse()`/`.as_matrix()` are
  FACTORIES returning DISTINCT objects (`_AdjointOperator` :704, `SweepOperator`
  `sweep_operator.py:44`, an ndarray). `A` and `A⁻¹` and `A.H` are three
  DIFFERENT morphisms sharing the substrate — NOT "views of one operator."

VERDICT on the user's hypothesis:
- "loss representation is the shared substrate, NOT an operator" — CONFIRMED,
  and ALREADY REALIZED (`InvertibleOperator` `streaming.py:520` routes apply
  :732 / apply_transpose :760 / solve :797 all into the ONE `loss_representation`
  property :703; `.inverse()` :781 returns `SweepOperator(self)`).
- "distinct OPERATOR VIEWS each a LinearOperator with ONE promise" — CONFIRMED
  (this is the RIGHT model).
- Q1 SLOGAN "forward-substitution and forward-matvec are different VIEWS of ONE
  operator" — **FLAWED, broken.** `A` and `A⁻¹` are TWO operators sharing ONE
  substrate; the SUBSTRATE has the views, the operator pair does not. This exact
  slogan is what licenses the `_GaussSeidelResolvent` confusion (apply+solve on
  one object). Correct phrasing: "two operators, one substrate," never "two
  views, one operator."
- Q5 "`as_matrix()` is a co-equal FOURTH view" — **PARTIALLY FLAWED.** Forward/
  Adjoint/Inverse are ENDOfunctors (return `LinearOperator`); `as_matrix()` is a
  functor OUT of the category (returns ndarray, the serialization boundary,
  dimension-gated). NOT co-equal — a different KIND of arrow. The co-equal fourth
  is the dense INVERSE strategy `DenseInverseOperator`, and it already exists
  un-named as `_as_dense`+`direct_eigenvalue` (`homogeneous/solver.py:119/129/195`).

## The four inverse STRATEGIES (user under-counted at 3)

`A.inverse()` is a STRUCTURE-KEYED FACTORY (+ optional `strategy=` override):
- triangular `(L+C)` → `SweepOperator` (direct fwd-substitution; LU frame)
- general `(L+C−S)` → `GreenOperator` = `(L+C−S)⁻¹` = Neumann series of the
  (L+C)-preconditioned splitting. **SourceIteration IS this via Richardson;
  KrylovAccel via GMRES** — so GreenOperator is the OPERATOR FACE of the SI/Krylov
  loop and MUST WRAP the driver (Smell #16 shape-4 if re-implemented). DEFERRED
  (no literal consumer yet — k-path is power-iter-with-SI-inner). The L-007
  resolvent backbone `fix(step)` predicts this exactly.
- small / CP `[P]` / 0-D → `DenseInverseOperator` (`lu_solve([A],·)`); already
  IS `direct_eigenvalue`'s dense `A⁻¹F`, un-promoted.
- parametrized `(sT+A)` → `ResolventOperator(z)` (α-work, §22 shift-invert).

Adjoint asymmetry (Q2/Q6, predicted by triangularity): the DIRECT (sweep)
inverse's `.H` needs a NON-TRIVIAL reverse-DAG traversal (#280, `xfail`); the
KRYLOV inverse's `.H` is FREE (`GreenOperator(A.H)` = GMRES on transpose-matvec).
This is WHY the carve defers `SweepOperator.H` to #280 — a property of WHICH
algorithm, not an oversight. (Matches `streaming_apply_transpose_frame`.)

## Q4 verdict: subclasses, structure-keyed factory, NOT methods/free-strategies

`.inverse()` returns SUBCLASSES of `LinearOperator` (grand report §3.3 lists
Sweep/Green as direct siblings — no separate `InverseOperator` ABC needed),
selected by the operator's STRUCTURE (triangular⇒Sweep; general⇒Green;
small/CP⇒Dense), with a strategy-override SEAM for the exceptions (force Krylov
on triangular for DSA; force Dense on a block). NOT (a) methods on one class
(adjoint structure differs per kind — would replicate the resolvent confusion);
NOT pure (c) free strategies (a CP matrix can't sweep, a meshed (L+C) can't
densify — the operator OWNS the selection).

## NEW SMELL (Part C promotion candidate): "layer-confusion / vestigial-forward"

An OPERATOR-layer object (a LinearOperator = ONE morphism, ONE promise) carrying
SUBSTRATE-layer multiplexing — an `apply` AND a `solve` that are inverses of
DIFFERENT operators. THE WITNESS: `_GaussSeidelResolvent` (`sn/solver.py:338`)
is duck-typed (NOT a LinearOperator), `apply` :386 = `(L+C)ψ` documented "NEVER
called", `solve` :391 = `(L+C−B_lower)⁻¹q`. The vestigial forward exists ONLY to
satisfy a duck-typed `CAP_APPLY` slot. TELL: an `apply` a docstring admits is
never called, on an object whose `solve` inverts a different operator. FIX: it
is a `SweepOperator` on the boundary-folded `(L+C−B)`, carrying NO `apply` of its
own (an inverse operator's `apply` IS the inverse action); the vestigial forward
deletes when `SourceIteration` consumes `L.inverse().apply` (carve P3).
DISCRIMINATING test = round-trip `inv.apply(A.apply(x))==x`: `_GaussSeidelResolvent`
gives `(L+C−B)⁻¹(L+C)x ≠ x` (off by the B reflection) ⇒ REDs; a clean
`(L+C−B).inverse()` passes. The round-trip is the structural PROOF the type is
confused. Distinct from Smell #16: that is two REPS of one quantity; this is one
operator-object faking a forward it never uses to inhabit the wrong layer.
PROMOTION STATUS: first sighting. Held for a second (a non-resolvent
operator-layer object carrying an unused method to satisfy a duck-typed slot)
before skill Part C; fire inline until then.

## The four foreign frames that fired (all concrete, all with fail-able tests)
- **Dagger category w/ partial inverses** (inverse/restriction category): `.H`
  involutive endofunctor, `.inverse()` partial (only isos), `is_invertible` =
  the restriction idempotent; `(A.H).inverse()==(A.inverse()).H` is FUNCTORIALITY
  (a theorem = #280's witness). The frozenset→predicate carve IS the idempotent
  made explicit.
- **LU / triangular forward-substitution**: `(L+C)` is lower-block-triangular in
  DAG order; matvec and inverse are the SAME factor traversed opposite ways;
  adjoint-inverse = reverse-order substitution on the transposed factor.
- **Resolvent formalism** `(A−z)⁻¹`: GreenOperator = Neumann series of the
  preconditioned splitting; SI/Krylov = its application strategies; `(sT+A)⁻¹`
  = the parametrized α-family.
- **Matrix-free↔explicit duality**: `as_matrix()` = materialization functor
  `Op→Mat` (apply-to-basis, `_as_dense`), dimension-gated, a functor OUT.

## Cross-method pollination (durable)
- MoC: ray-trace solve is ALSO triangular fwd-substitution (along
  characteristics) ⇒ `TrackSweepOperator`, sibling of `SweepOperator` — the
  SECOND witness that earns SweepOperator a generic `numerics` home (user already
  flagged the relocation trigger).
- CP: `[P]` is DENSE by construction ⇒ CP lives natively in the
  `DenseInverseOperator`/`as_matrix()` branch — CP is the consumer that PROMOTES
  `as_matrix()` from a 0-D homogeneous convenience to production machinery.
- Diffusion: self-adjoint elliptic, NO triangular factor ⇒ `.inverse()` is a
  CG/multigrid GreenOperator and `.H==self` (the L-007 diffusion EXCEPTION) —
  the negative control proving the factory keying isn't SN-overfit (exercises
  Green+Dense, never Sweep).
