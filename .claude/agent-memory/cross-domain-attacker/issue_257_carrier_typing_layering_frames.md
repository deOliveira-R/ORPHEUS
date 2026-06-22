---
name: issue-257-carrier-typing-layering-frames
description: Frame attack on the #257 carrier-type-hierarchy + numerics↛transport layering design — confirms Vector Protocol is irreducible (forgetful-functor image), FullField/TimedFullField split forced (Cofree comonad), Field-relocation orthogonal hygiene.
metadata:
  type: project
---

# #257 carrier-typing + layering design — 6-frame attack (CONFIRM explorer)

DESIGN MEMO (not bug). Question: is the `numerics ↛ transport` (L1↛L2)
barrier — hence the structural `Vector` Protocol apparatus — an artifact
of a bad folder hierarchy (user hypothesis), or irreducible (explorer)?
Plus: rename mis-named `TransportState` Protocol; #217 split timeless
`FullField` out of `TimedFullField`. VERDICT = option (A): keep layering,
rename to timeless `FullField`, history orthogonal. CONFIRMS explorer and
upgrades "irreducible" from assertion to theorem.

## Ground truth verified (not just prose)
- `Field` (`orpheus/numerics/field.py`) has EXACTLY ONE caller:
  `test_field_is_abstract_direct_instantiation_rejected`. Nothing in
  numerics constructs/consumes it → explorer's relocation claim CONFIRMED.
- `numerics` references `transport` ONLY in docstring `:class:`
  cross-refs (vector.py, iteration.py, full_field_space.py), NEVER in
  imports. iteration.py:170/240 explicitly: "numerics MUST NOT import
  transport."
- Barrier is `FORBIDDEN_EDGES["numerics"] = L2|L3` in
  `tests/test_layer_imports.py:55` — an enforced parametrized test.
- Biproduct ALREADY named + implemented: `operator.py:109-120` types
  operators as 2×2 block matrices "by the biproduct theorem";
  `FullFieldSpace = V_bulk ⊕ V_trace` with block-diagonal metric
  `G = G_bulk ⊕ G_trace` (`spaces/full_field_space.py`).
- `Vector` Protocol (`numerics/vector.py`): 4 dunders
  `{__add__,__sub__,__rmul__,__truediv__}`, `V=TypeVar(bound=Vector)`,
  satisfied by ancestor-disjoint `np.ndarray` / `Field` leaves /
  `TimedFullField`. Docstring "documented falsehood" = the type was
  widened from the `np.ndarray` lie to the structural truth.

## The 4 confirming frames (CONCRETE PAYOFF — promotion candidates)

### F1 — Forgetful functor / free-forgetful adjunction (STRONGEST)
`Vector` IS the object-image of the forgetful functor `U: C_carrier →
C_vec` (forgets units/mesh/space-identity, keeps the 4 vector ops).
`V=TypeVar(bound=Vector)` = "natural in the U-image of any carrier."
The barrier = `U` has no nameable inverse from inside its codomain
(`numerics` lives in cod(U)). PROMOTES "irreducible" to a THEOREM: the
Protocol is the SHADOW of the layer DAG, not an independent choice. Only
lever to dissolve `Vector` = move the algebra UP to L2 (abandon
method-agnosticism). Field-relocation can't touch `Vector` — Field is
not on U's path. First test: build a carrier satisfying `Vector` that is
NOT Field/TimedFullField/ndarray (moment-only iterate, tuple-of-fields);
if algebra+drivers consume it unchanged, U-image confirmed.
TRIGGER (new candidate for A.2): "abstract algebra at layer N acting on
carriers defined at layer N+1, bridged by a structural contract" →
forgetful functor; the structural Protocol IS U's image, irreducible by
the layer DAG.

### F2 — Fibration / parametric polymorphism over carrier category
`p: E(transport carriers) → B(vector spaces)`; `apply: V→V` is a
CARTESIAN morphism (base-defined, lifts uniquely per fibre because
operators are role-preserving — flux→source role-change is INSIDE the
fibre, vector-shape fixed). The "layer tension" the user fears IS the
defining property of a fibration (base ↛ total space). CHALLENGES the
QUESTION's framing (barrier is not an artifact), CONFIRMS explorer.
Sharpens: single-`V` (vs rejected per-op `Generic[F]`, #226) is correct
BECAUSE operators are fibre-preserving — load-bearing, not cosmetic.
⚠ FIRST TEST = the one place a frame pushes back on CODE: does any
operator `apply` branch on `isinstance(x,...)`/`type(x) is` rather than a
structural `.bulk`/`.boundary` accessor? If yes, `Generic[V]` is a
partial lie + fibration leaks. Block-role dispatch is the suspect. CHEAP,
do before #217 hardens carrier types.

### F3 — Module over a ring / multiplier algebra (#257 already names C=M[σ_t])
Carrier ψ = module M over ring R (coefficient fields); operators =
`End_ℝ(M)`. Classical fact: `End_R(M)` needs only M's additive group +
R-action, NOT M's concrete realization → endomorphism algebra rightly
at the `Vector` level (additive group + scalar action). RELOCATES the
user's hypothesis: the interesting filing question is about R (the
coefficient RING), not M (the carrier). `Field` is OVERLOADED — additive-
group base of flux-module AND (via CrossSectionField) underlying set of
the scalar ring. Field-relocation + #257 coefficient-algebra split are
the SAME hygiene act from two angles — DO TOGETHER. (Aligns with
[[coefficient-field-promotion-frames]] F2: CoefficientField is a
commutative algebra, not a flux role.) First test: do composers need the
RING at construction (typed CrossSectionField) or only raw arrays at
apply? `DiagonalOperator(weights:ndarray)` is ring-agnostic → End below
ring → L1 correct.

### F4 — Cofree comonad (the FullField vs TimedFullField split, #217)
`TimedFullField = Cofree(FullField, depth=d) = FullField × FullField^d`.
`extract = at_lag(0)`, `advance (timed_full_field.py:329) = comonadic
extend`, rotating buffer = the cofree tail. A comonad W acts on the
category of TIMELESS objects; operators are base arrows FullField→
FullField, lifted through W ONLY by iteration drivers. The operator
algebra NEVER needs W. FORCES #217: operator output (`Cψ` source, `b−Ax`
residual) is a base object — typing it `TimedFullField` gives a vestigial
history tail = type error of ALTITUDE (Smell #16 shape-3 at one remove:
output typed with the iterating-state's decoration, like Δψ-typed-as-
state). First test: zero production sites read `.at_lag(k>0)` on an
operator RETURN value (all history consumers are drivers). NAMING: buffer
serves time AND Krylov ("time/Krylov stencil" — vector.py:14,
iteration.py:184) → "Timed" under-describes → rename
`TimedFullField → IteratedFullField` (or `FullFieldHistory`); if kept,
document as `Cofree(FullField,d)`.

## Frames that fired but were already-matched / cosmetic-resolving
- Biproduct (A.2): already named+implemented in code. PRESCRIBES for the
  open Qs: `FullField` is ONE object (biproduct = product AND coproduct,
  `ι_b π_b + ι_s π_s = id`), not "a pair stored together"; timeless
  `FullField` should BE the biproduct, history = comonad ON TOP. New test
  it suggests: round-trip split-then-rejoin == id (catches bulk/boundary
  overlap/gap, e.g. double-counted tangential ordinate — a real bug).
- Structural-vs-nominal subtyping (A.2): `Vector` structural = the
  bounded-quantifier UPPER BOUND (`∀V<:Vector`); nominal `FullField` =
  the concrete instantiation. They COEXIST, different syntactic roles.
  Bound MUST be structural (admits ndarray, no shared ancestor); carrier
  base CAN be nominal (transport leaves are all ours). The user's worry
  "Protocol is a hack for a missing base class" is itself the error — the
  Protocol is the quantifier bound, not a stand-in for FullField.

## Cross-method pollination
- Coupled-multiphysics: `(ScalarFlux, TemperatureField)` product carrier
  should satisfy `Vector` (componentwise dunders). First test: drive a
  2-component product through `power_iteration` unchanged → validates
  `Vector` as universal carrier interface (the campaign's central claim).
  If it needs `.copy()`/`.to_flat()`, Protocol under-specified — close
  the gap BEFORE multiphysics coupling lands.
- Krylov-vs-time stencil: same history-window structure (vector.py says
  "time/Krylov stencil"). F4 predicts they unify (Krylov basis {ψ,Aψ,…}
  IS a history window). First test: do both paths consume
  `.history` through one accessor or two? Two = Smell #16 shape-1
  (one history concept, two reps) → fold onto one `Cofree(FullField,d)`.
  Confirms the F4 rename (buffer is method-neutral).

## UNEXPLORED (no trigger)
topology/manifold (no geometry-variation), group theory (no carrier
symmetry), tensor-network/MPO (biproduct rank-2 FIXED not a bond dim —
distinct from Variant-α MPO match), QMC (no integrand), homogenization
(no scale separation), sheaf/gluing (blocks DISJOINT not a cover —
biproduct already captures it, sheaf adds nothing).

## Memory-discipline note
F1 (forgetful-functor → structural-Protocol-is-irreducible) and F4
(Cofree comonad → timeless-base-split forced + Smell-16-shape-3-at-one-
remove for history-decorated operator output) are the strong promotion
candidates for skill Part A.2 at next revision. F2's `isinstance` leak
test is the one CODE action to carry to plan-mode.
