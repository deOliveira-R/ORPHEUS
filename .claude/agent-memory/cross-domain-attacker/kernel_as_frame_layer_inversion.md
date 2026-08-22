---
name: kernel-as-frame-layer-inversion
description: The Frame precedent read at HEAD — Frame is the 2-field BINDER, Basis is the rich DATA class; a "kernel should generate its operators" proposal is off by one layer, and the data object's verbs must return ARRAYS while only the binder returns OPERATORS.
metadata:
  type: project
---

# Kernel-as-Frame — the layer inversion, and the rule that decides it

**Fact.** `orpheus/numerics/frame.py`'s `FrameBase` is a `@dataclass(frozen=True)`
with **exactly two fields** (`basis`, `measure`) that implements **zero math** —
both operator faces' `apply` bodies are one-line delegations back to the basis
(`frame.py:434-437, 475-476`). The rich, math-bearing class in that stack is
**`Basis`** (`numerics/basis/base.py:117`): 6 representation-free verbs whose
representation-dependent arguments (`table`, `weights`) are supplied by the
binder at call time, plus `gram_structure` and `space`.

**Why:** a proposal of the form *"class X is thin; let it assemble its operators
the way Frame does"* maps X onto the wrong layer. A thin DATA class's analogue is
`Basis`; the assembler's analogue is the BINDER. Read that way the Frame
precedent **argues FOR an external binder** — `FrameBase(basis, measure)` IS
`bind(kernel, space)` with the nouns changed — while the *thinness diagnosis*
survives, re-aimed at the data class (`ScatteringKernel` has 1 representation-free
morphism; `Basis` has 6).

**How to apply:** on any "should this object generate its operators" dispatch,
map the problem's objects onto **both** layers of the cited precedent and keep the
mapping that preserves ARITY (field count, verb count). Then apply the
discriminator:

> ⭐ **In a data/binder split, the DATA object's verbs return ARRAYS; only the
> BINDER returns OPERATORS.**

That one rule is grep-checkable and it decides the layering question outright:
`basis.analyze(...) -> NDArray` keeps `base.py` at **zero runtime imports**
(stdlib + numpy; `DiscreteMeasure`/`FunctionSpace` are `TYPE_CHECKING`-only), which
is the single property making it reusable from any layer. A `datum.bind(space) ->
LinearOperator` inverts exactly that.

## The four durable findings from the CS4b stress test (2026-08-21)

1. **Binding is a BINARY operation, so neither operand owns it.** The precedent's
   answer is a **THIRD OBJECT** — not `basis.bind(measure)`, not the reverse. And
   the third object is where the caching lives, which neither operand can host
   *because neither knows the other*. ⟹ the chartered `EE-6` binding base IS
   `FrameBase`, and it is **one field short**: lift the DATUM into it
   (`BoundOperator(datum, space)`), and the abstract `data_ng` hook deletes.
2. **A "declare which analysis verb" requirement belongs to the FRAME TYPE, and
   this repo already ruled it** — `numerics/projection.py:24-36` retired
   `GalerkinProjection`/`PetrovGalerkinProjection` marker ABCs at #268 for
   declaring a discipline on the *role* when it belongs to the *frame*. Putting it
   on a datum is the same category error one layer out, and it is structurally
   impossible: `analysis` vs `project` differ by `G⁻¹`, and `G` is a row-sum probe
   over basis AND measure — a datum has neither.
3. **A 3-of-4 uniformity gap is INFORMATIVE, not a smell, when the fourth member
   is the degenerate case.** The collision multiplier `C` has no frame because a
   diagonal operator's eigenbasis is the nodal basis, whose frame is `Id`
   (L-009). And the asymmetry is already a TYPE with a gate:
   `IntegralKernelOperator`'s sole discriminator is the `kernel` member, so
   "uniformising" C would blind a gate that currently discriminates. Criterion to
   use: *different construction ⟹ smell; degenerate case of the same construction
   ⟹ informative, and unifying DELETES content.*
4. **`bind` is a DAGGER FUNCTOR, natural in the space — and that unifies three
   separately-chartered gates.** multiplicativity `bind(K₁K₂)=bind(K₁)bind(K₂)` =
   functoriality; `bind(K†)=bind(K)†` = the dagger law; `bind(condense K)=T·bind(K)·T⁺`
   = naturality. It also refutes datum-owned binding independently: **a functor acts
   on MORPHISMS too**, and `condense`'s image is an operator-level conjugation
   involving two spaces, which no single datum instance can own. (This is the rare
   case where the categorical frame earns its place past L-001's standing
   suspicion — it produces a refutation, not a restatement.)

## The unowned-fusion refutation (the other independent kill)

`ScatteringOperator.full_scatter_kernel` is `frame.conjugate(Λ_{ℓ≥0} + N2N)`
(`scattering.py:678-683`) — a binding of a **channel SUM**, taken in moment space
BEFORE conjugation because the project pins the scattering kernel at **0 ULP**.
`K_scatter.bind(space) + K_n2n.bind(space)` is a *different number*. ⟹ any
per-datum `bind` API has no owner for the fused operator. The reshape that gives
it one: the binder consumes an **emission kernel** `K_emit = Σ_c ν_c Σ_{r,c}`,
formed by data-level addition.

## ⚠ The god object is not hypothetical

`ScatteringOperator` at HEAD **IS** the "rich class generates its operators"
design: one class holding the datum + the representation and emitting 4 operator
factories, its own frame, 9 data read-throughs and 2 assembly helpers. The
mechanism to carry forward is not "big classes are bad" — it is
**the object that owns the generation set accretes the generation set.**
