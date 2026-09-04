---
name: transfer-family-carve-rulings
description: #426 step 2 (one transfer family, the (n,2n) gain anisotropic) — the four durable review rulings, the two new recurring smells (the role-body twin; the value-vs-shape predicate promotion), and the probes that settled them
metadata:
  type: project
---

#426 step 2 `1a3b78ec` (`feature/n2n-transfer-family`): two collision-gain terms
collapsed into ONE Legendre-stack-with-a-yield kernel (`TransferKernel`), one field,
one angular core (`TransferOperator`), one P0 core (`IsotropicTransfer`), with
`ScatteringOperator` / `N2NOperator` / `IsotropicScattering` / `IsotropicN2N` as thin
ROLE subclasses gated by an AST test. Reviewed 2026-09-04. Verdict: VIOLATIONS
REQUIRE REWORK (3, all prose-or-plumbing; the numerics were clean).

**Why:** the rulings below are about SHAPES that recur across ORPHEUS carves, not
about (n,2n). **How to apply:** read before any review of a core/role split, a
value-predicate promotion, or a "two terms of one algebra" collapse.

## The four durable rulings

**1. A role's content is DATA, not methods — and a `ClassVar` is how you say so.**
The carve got half of it right (`isotropic_binding: ClassVar[type[...]]` naming the
role's P0 binding) and half wrong (four hand-written extraction classmethods
differing in ONE token: `TransferMaterialField.scattering` vs `.n2n`). The repair
that generalises: put a `channel: ClassVar[Callable[[Facade], Field]]` on the role
and ONE mint on the core. `[M]` verified — a bound classmethod stored as a `ClassVar`
is NOT re-bound, is NOT a dataclass field (Py 3.14.3), and reads identically off
class and instance. ⭐ The bug-habitat argument that makes this a VIOLATION rather
than a nit: **the defect the carve retired WAS the two mint bodies having diverged**
(N2N minted at `L=0` while S minted at `scattering_order`), and the carve repaired it
by editing BOTH bodies — i.e. it fixed the symptom and re-shipped the habitat. When
you see a carve retire a twin, ask whether the twin's *generating mechanism* went
with it.

**2. A shape predicate promoted to a VALUE predicate is honest iff the two branches
are measurably the same VALUE — and then the old spelling becomes a stale synonym.**
`scattering_order == 0` → `is_isotropic` (every moment above ℓ=0 exactly zero) was
forced: once both bindings run at the solve's order, "order 0" stops naming the
isotropic case. `[M]` the skip is a PURE short-circuit — `array_equal` True,
`max|Δ| = 0.0`, 0 signed-zero differences, measured by forcing the non-skip branch
in-process. ⚠ The residue to hunt every time: the OLD spelling survives elsewhere and
now means something different. Here `kernel`'s guard and the crosscheck's fixture
assertion both still say `scattering_order == 0/≥1` while CALLING it
"isotropic"/"anisotropic". A promotion creates synonyms; grep the word, not the
symbol.

**3. Deriving a sub-binding beats storing it — and the guard you don't need is the
proof.** `TransferOperator.isotropic_energy` DERIVES the P0 binding from its own
datum at order 0 (a `cached_property` + a role `ClassVar`), so the two bindings of
one channel cannot disagree. `FissionOperator.energy` STORES it as a field and pays
with a runtime shape guard (`fission.py:276-285`). Same algebra, two shapes, and the
derived one is strictly stronger. ⟹ when reviewing "should X be a field or derived",
the tell is whether the alternative needs a validator.

**4. `at_order` returning `self` at the identity + routing both other arms through
`dataclasses.replace` is the textbook Pattern 4 ∩ 2** — the invariant lives once in
`__post_init__` and every same-type-producing path re-runs it. `[M]` verified,
including that padded buffers are distinct and read-only (no shared zero alias).
Copy this shape.

## Two smells this carve added to the catalogue

**S-A. A guard whose COMMENT names two things and whose BODY checks one.**
`TransferOperator.__post_init__` says *"both minted faces must be bound to THIS
binding's interior"* and checks `flux_analysis` only. `[M]` measured live: a
`source_reconstruction` from an 8-ordinate quadrature on a 4-ordinate interior is
ACCEPTED, and the windowed arm then returns an `(8,2,4)` source **with no raise**
(the scalar `iso/W` broadcasts). Pre-existing, carried verbatim into the new shared
core — so the exposure doubled while the text stayed. ⟹ on any `__post_init__`,
count the nouns in the comment and the clauses in the body.

**S-B. A gate whose DOCSTRING cites the measurement that its BODY ignores.**
`test_transfer_roles.py`'s population row exists to answer "is this the whole
population?", its docstring quotes the direct-vs-recursive subclass measurement, and
its body calls `__subclasses__()` — direct only. A subclass of a ROLE is invisible to
both filters. ⟹ when a gate's docstring supplies a number, check the body uses it.

## Probes that settled it (reusable; all in-process, ~seconds)

* **value-vs-shape short-circuit**: build the operator twice, monkeypatching the
  value predicate to `False` on the FIELD class so the full product runs on zeros;
  compare with `array_equal` AND `np.signbit` (a `-0.0` difference is the one thing
  `array_equal` reports as equal but a bit-comparison would not).
* **admission census**: loop the illegal values `(0, -1, True, False, 2.0, 1.0, "2")`
  — `True`/`False` is the load-bearing pair (`bool ⊂ int`).
* **is-it-a-dataclass-field**: `{f.name for f in dataclasses.fields(C)}` settles the
  ClassVar-on-a-dataclass question in one line. (Confirms the CS4a refutation: under
  `from __future__ import annotations`, Py 3.14 DOES detect a stringized `ClassVar`
  on a non-frozen dataclass. Do not re-derive this from AGENT.md's older note.)
* **mismatched-face admission**: build the operator with a face from a second mesh's
  frame and APPLY it — the shape of the output is the finding, not the construction.
