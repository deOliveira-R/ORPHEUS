---
name: issue-256-vector-protocol-step1
description: #256 Step 1 Vector Protocol — CHANGES-REQUESTED; the "honest minimal dunder set" audit MUST include the eigenvalue.power_iteration driver, not just operator.py + iteration.py
metadata:
  type: project
---

# #256 field-typed operator algebra — Step 1 `Vector` Protocol review

CHANGES-REQUESTED (branch `feature/field-typed-operator-algebra`, uncommitted). Step 1
mints `orpheus/numerics/vector.py`: a `runtime_checkable` structural `Vector` Protocol
(`__add__`/`__sub__`/`__rmul__`) + `V = TypeVar("V", bound=Vector)` so a LATER step can
retype the algebra's `apply(self, x: V) -> V` (the honest endomorphism). No consumers
change in step 1 — it only DEFINES the contract. Conforming carriers: `np.ndarray`, every
`Field` leaf, `TimedFullField`.

**Why:** the campaign's stated acceptance criterion is "the dunder set is the honest MINIMAL
set the operator algebra actually invokes." Step 3's `apply(x: V) -> V` retype only
type-checks if `V` carries every op the bodies use.

**How to apply (the load-bearing lesson):** when auditing "what does the algebra invoke on a
*vector*", the driver surface is operator.py + iteration.py + **eigenvalue.py**. `KEigenvalue`
(iteration.py) DELEGATES its outer loop to `power_iteration` (eigenvalue.py:1119-1120), so
eigenvalue.py IS part of the same iteration-driver surface. The author scoped the audit to
operator.py + iteration.py and MISSED two `vector / scalar` (`__truediv__`) sites:
- `iteration.py:1021` `self.F.apply(flux_distribution) / keff` (k-posing eigen-source Fψ/k)
- `eigenvalue.py:226` `flux_distribution / p` (per-step production renorm, ERR-052)
Both carriers (`Field.__truediv__` field.py:279, `TimedFullField.__truediv__`
timed_full_field.py:319) HAPPEN to provide it, so the gate doesn't catch the omission — but
the Protocol is then a "documented falsehood" of exactly the kind vector.py:13 indicts, and
the Step-3 retype of `Fψ/k` will fail to type-check (or get silenced with an ignore = the
falsehood ships). BLOCKING fix: add `__truediv__(self, scalar: float) -> Vector` + a
`test_no_scalar_div_rejected` mirror.

**What was RIGHT (reinforce):**
- The `__rmul__`-over-`__mul__` correctness fix is exemplary: scalar-mul in the composers is
  exclusively LEFT-mul (`scalar * v`) — `ScaledOperator.apply` operator.py:867,
  `ZeroOperator.apply` operator.py:955 (`0.0 * x`). Docstring cites both sites (Pattern 3).
  The additive identity `0.0 * x` (no abstract `zero` member needed) is a genuine elegance win.
- Protocol-vs-ABC-vs-Generic[F]: structural Protocol is CORRECT. np.ndarray can't inherit an
  ORPHEUS base (kills ABC); the inner leaf changes role (flux→source) at runtime while the
  CARRIER type is fixed (kills per-operator Generic[F]). Mirrors the existing `LinearOperator`
  / `EigenvalueSolver` Protocol idiom. Promotes the ad-hoc `_is_ravellable` duck-type
  (iteration.py:178) to a named contract — right direction.
- Axis-4 NON-issue: Protocol members returning the wider `Vector` (not `V`) is FINE for
  `apply(x: V) -> V`. The concrete narrowing flows through the carriers' self-typed dunders
  (`Field.__add__` `self: T -> T` field.py:262); `Vector` is the correct upper bound.

**STANDING TELL — flat-axis-primitives are NOT carrier-Protocol consumers.** `DiagonalOperator`
/ `RankOneOperator` apply bodies do `self._reshape() * np.asarray(x)` — they coerce to ndarray
and act on ONE tagged axis of `.values`, NOT on the typed carrier through Vector dunders
(vector.py:44-46 already documents this). So they impose NO `Vector` requirement. Don't
mistake an operator that internally touches `.values` for a carrier-level vector op.

**Out-of-scope (correctly): `.copy()`** — invoked on the carrier in `power_iteration:203`
(`flux_old = flux_distribution.copy()`) but it's an aliasing utility, NOT a vector-space op,
and never appears in an `apply(x: V) -> V` body. Recommend a docstring sentence scoping
`Vector` to vector-space ops only (`.copy()` / `to_flat`/`from_flat` stay on `_is_ravellable`)
so a future reader doesn't "complete" the Protocol with `.copy()` and narrow what conforms.
