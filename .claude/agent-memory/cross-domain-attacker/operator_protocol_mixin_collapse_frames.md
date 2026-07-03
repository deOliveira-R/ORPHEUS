---
name: operator-protocol-mixin-collapse-frames
description: Frame attack on the P4.5a two-param LinearOperator(Protocol[Din,Cout]) / LinearOperatorMixin(Generic[V,W]) split. VERDICT the stateless dunder-installer Mixin is Smell #16 shape-1 (two co-required contracts for ONE concept, bridged by a TYPE_CHECKING apply stub at operator.py:446 that exists ONLY to reconcile the two halves). Two-param change ENABLES the collapse. Frame-1 = dunders-as-default-methods ON the Protocol (operators are morphisms; the hom-set carries the algebra) — GATED on a pyright spike (can a contravariant Din host construction-returning dunders __add__->OperatorSum[Din,Cout]?). Frame-2 fallback = single Operator(Generic[Din,Cout],ABC), drop variance (the composers ALREADY use invariant V/W — declared Protocol variance is unused below the read surface; operators are never np.ndarray so the structural-Protocol justification is miscopied from Vector). Frame-3 (fibered category) RULES the leaf-vs-composite question: type the LEAF [Flux,SourceSink] (genuine cartesian morphism), composite [FullField,FullField] is DERIVED from source-fiber closure, NOT a primary contract. RETIRE-list ranked: 7 @overload stubs (collapse to one generic apply, singledispatch already does it) > Mixin apply stub > Mixin class > variant TypeVar machinery + _AdjointOperator [W,V] workaround. Branch-verified refactor/operator-inverse-algebra 2026-06-24.
metadata:
  type: project
---

# Operator Protocol/Mixin collapse — structural frame attack (P4.5)

Branch-verified DIRECTLY on the worktree (L-005), NOT Nexus. Files: `orpheus/numerics/operator.py`
(Protocol `:274`, Mixin `:371`, the `TYPE_CHECKING` stub `:436-446`, composers), `orpheus/numerics/vector.py`
(Vector/V), `orpheus/sn/scattering.py:1388-1408` + `orpheus/sn/fission.py:466-482` (the `@overload` blocks).
This is the typing-LAYER lift of [[flux_to_sourcesink_operator_contract_frames]] (which ruled the
CONTRACT = linear part of an affine map Flux→SourceSink): that memo said "type the Operator[Domain,Codomain]
contract"; THIS memo rules HOW the base class should carry it under PEP-696 two-param.

## The ground-truth structure (the defect)
- `LinearOperator(Protocol[Din,Cout])` — the READ contract, VARIANT (Din contra, Cout co). `:274`.
- `LinearOperatorMixin(Generic[V,W])` — the dunder-installer, INVARIANT V/W. STATELESS (no apply, no
  __init__). `:371`. Carries a `TYPE_CHECKING`-only `apply(self,x:V)->W` stub `:436-446` whose IN-CODE
  comment says it exists ONLY "to let `self` satisfy the LinearOperator Protocol inside the algebra dunders."
- A concrete operator must satisfy BOTH: inherit Mixin (for `+−*@.H`) AND structurally match Protocol
  (for `isinstance`/`_has`). TWO inheritance obligations for ONE concept = Smell #16 shape-1.
- Variance declared on the Protocol does NOT propagate: 15/18 concrete leaves inherit `LinearOperatorMixin`
  BARE (no subscript); only `StreamingOperator["FullField"]`, `MultiplicationOperator["FullField"]` + the 6
  composers subscript it. The honest `Din≠Cout` is patched per-leaf by `@overload` (Smell #16 shape-2:
  signature hand-copied in two places, dropped to bare in between).
- `ScatteringOperator`/`FissionOperator` inherit bare, then `@overload` `apply` in a `TYPE_CHECKING` block
  CONFESSING "is NOT an endomorphism V→V (the mixin's nominal contract)". Runtime `apply = _apply_impl`
  (`functools.singledispatchmethod`). 4 scattering + 3 fission stubs.

## The three native frames (each FAIL-ABLE)
- **Frame 1 — operators are morphisms; put the dunders ON the Protocol as default methods.** Python Protocols
  MAY carry concrete method bodies; a host inherits them. Then the `:446` stub DELETES (self IS a
  LinearOperator by construction) and the Mixin retires. DISCRIMINATOR (REDs on current bare-Mixin):
  `class Probe(LinearOperator[Flux,SourceSink])` WITHOUT inheriting the Mixin — `(Probe()+Probe())` must
  return an OperatorSum; today RAISES AttributeError __add__ (dunders only on Mixin). ⚠ GATE: a
  construction-returning dunder `__add__->OperatorSum[Din,Cout]` puts the contravariant Din in
  invariant-construction position — the SAME wall P4.5a hit (operator.py:76-83). NEEDS A PYRIGHT SPIKE
  before committing. Spike outcome forks the plan: pass⇒Frame-1, fail⇒Frame-2.
- **Frame 2 — single `Operator(Generic[Din,Cout],ABC)`, drop variance.** `apply`=`@abstractmethod`, dunders
  concrete on the same class, domain/codomain/block_role/capabilities concrete defaults. Variance dropped:
  the composers ALREADY use invariant V/W TODAY, so declared Protocol variance is UNUSED below the read
  surface. KEY structural attack: the structural-Protocol justification (ndarray conforms without
  inheritance, vector.py:45-68) is CORRECT for Vector (the carrier) and MISCOPIED onto LinearOperator (the
  operator) — operators are NEVER ndarray, every operator is an ORPHEUS class already inheriting the Mixin,
  so the operator's Protocol-ness is unmotivated and the Mixin is the patch the no-default-methods stance
  forced. DISCRIMINATOR: a subclass omitting `apply` RAISES at class-definition time (abstractmethod), not
  call time — moves the failure earlier; current bare-Mixin cannot.
- **Frame 3 — fibered category; type the LEAF, derive the composite.** Carrier fibered over "physical
  quantity"; flux/source role = fiber coordinate; `S:AngularFlux→AngularSourceSink` = cartesian morphism
  (moves between fibers over fixed base = same mesh/group). Composite `(L+C−S−F):FullField→FullField` is
  endomorphic ONLY as a DERIVED property of source-fiber closure (sources sum closed) — NOT a primary
  contract. RULES the brief's leaf-vs-composite question: leaf role-change is PRIMARY, composite
  endomorphism is DERIVED. Typing the composite `[FullField,FullField]` as primary INVERTS the dependency
  and is WHY the role-changing leaves need `@overload` patches. DISCRIMINATOR: `(S+F).codomain` is
  AngularSourceSink and `(S+F).apply(ψ)+(S+F).apply(ψ)` type-checks while `ψ+(S+F).apply(ψ)` RAISES; a
  composite-primary design lets the cross-role sum through (both "FullField").

## RETIRE-list (ranked confidence×payoff)
1. The 7 `@overload` stubs → ONE generic `apply(x:Din)->Cout` (singledispatch already does the fiber
   dispatch; only the static surface is transcribed). HIGH.
2. The Mixin `TYPE_CHECKING` apply stub (`:436-446`). HIGH (deletes under Frame-1 OR Frame-2).
3. `LinearOperatorMixin` class itself. MED-HIGH, gated on the variance spike.
4. Variant Din/Cout TypeVars + 70-83 rationale + vector.py:78-83 + the `_AdjointOperator` `[W,V]` PEP-696
   ordering workaround (`:646-655`, "the ONLY composer with [W,V] order"). MED — the deliberate cost of
   dropping variance; retire only if the spike confirms variance unused.
5. (debt-payoff, not retirement) 15 bare-Mixin leaves gain honest `[Flux,SourceSink]` at the class header.

THE BLOCKING UNKNOWN = the pyright variance spike (Frame-1 vs Frame-2 fork). Run it FIRST.

## Cross-method (durable)
- Eigenvalue `fix(step)` layer ([[power_iteration_vs_keigenvalue_morphism]]): generality flows toward the
  OPAQUE interface — the operator base's public surface must be a SUBSET of what the iteration drivers call
  ({apply,solve,+,−,*,@,.H,capabilities,domain,codomain}); any base member outside is dead ceremony. Do NOT
  grow the base an `__init__` the drivers never touch.
- MoC fiber-bundle carriers (project memory "MoC structure"): the flux→source `1/cm` unit-gain (Frame-2 of
  the parent memo) = MoC's bundle transition function; type the units-gain ONCE in numerics so SN `S.apply`
  gain and MoC ray-segment transition read one place. Discriminator: gain readable off `S` WITHOUT applying
  (the parent memo's Frame-2 test).

## Refuted (durable UNEXPLORED)
- Homology — `.H` dagger pair is not a differential, no `∂²=0` (L-001).
- MPO/tensor-networks — TP/SumOfTP fixed-rank, no bond knob (L-001).
- Group/rep theory — SO(3) on quadrature/Λ, not the typing contract.
- Diff-geom/Christoffel — no curvature in the flat typing layer.
- Category-theory-as-abstract-nonsense — the categorical CONTENT (morphism/dagger/biproduct/fibration) IS
  the concrete frame and produces typed tests; no operad/PROP/2-cat lever adds a discriminating test (L-001).

Cross-refs: [[flux_to_sourcesink_operator_contract_frames]] (the CONTRACT verdict, parent),
[[issue_208_operator_algebra_frames]] (dagger biproduct category, the algebra spine),
[[issue_257_carrier_typing_layering_frames]] (apply:V→V is a fibration cartesian morphism — Frame-3 base),
[[power_iteration_vs_keigenvalue_morphism]] (generality→opaque interface, the cross-method rule),
[[issue_226_container_algebra_design]] (the structural Vector Protocol + apply(x:V)->V history).
