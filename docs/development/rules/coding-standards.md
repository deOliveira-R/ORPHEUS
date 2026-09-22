---
harness:
  kind: rule
  budget_tokens: 3300
  paths:
    - "orpheus/**"
    - "tests/**"
    - "derivations/**"
    - "tools/**"
  brief: >-
    tests run as `python -O -m pytest`; a bare `assert` outside a collected test module is stripped under `-O`, so a contract is a `raise` and a test-side check is `np.testing.assert_*`.
---

# Coding standards — the minimum-quality floor

Minimum standards every contributor, main agent and sub-agents, follows by default; Cardinal Rule 2 ranks them against `coding-elegance`. Founding cases: [the evidence page](../evidence/coding-standards.md).

## Clean before extending

Before adding a capability to a class/module, run a cleanup pass on that layer first: collapse double paths, move concepts to their native place, delete dead shims, fix twin sources of truth. The capability then lands as a **no-op extension through the one generic body**, not a third arm grafted onto debt.

- check: a plan proposing a capability extension inserts a **cleanup phase before** it; order the findings into those that must precede the extension, those that are independent polish, and those that explicitly wait; gate each cleanup substep bit-identical where possible.
- tell: the new arm needs a matching arm in the converter AND the constructor AND the gate. [case](../evidence/coding-standards.md#2026-06-11-from-axes-roundtrip)

## Type vs property — before minting a type

Mint a **type** iff (a) the concept has **two or more non-isomorphic realizations** AND (b) a **non-identity morphism** is actually applied to it. Otherwise it is a **property**: a field or flag on an existing type. One realization plus an identity change-of-basis makes a "type" theatrics: a conversion seam and no illegal state made unrepresentable (a single-basis spatial moment is a `property`, not a `SpatialOrder` type).

- check: count the realizations and name the morphism; identity-only change of basis means property.
- **An axis that changes the ARITHMETIC INTERFACE cannot be a phantom type parameter.** `Generic[Tag]` is erased at runtime and does not specialize dunders, so every instantiation shares ONE `__add__`; a torsor `A×V→A` that must forbid `A×A` and a vector `V×V→V` cannot share a body. Arithmetic or shape changes: a class; neither: a phantom parameter is allowed.
- tell: an implementation that "passes" only by branching on a stored tag at runtime is stringly-typed dispatch; `replace(obj, tag=Other)` type-checks and walks through the gate the type was minted to be.

## Naming — a name is the index that reconnects a fresh context to the right structure

A name is chosen for the signal it carries about the structure it names, and the agent judges it (`[R]` the user, 2026-06-23: "as long as it carries a strong signal"). Founding cases: [the evidence page](../evidence/coding-standards.md#2026-06-22-naming-rulings).

- check: name the parts in the object's native vocabulary, precisely: a discrete frame's faces are `analysis` (`M`) and `reconstruction` (`R`, the dual-frame synthesis), not "projection" and not `synthesis` (frame theory's naked synthesis operator, the adjoint of analysis over the SAME family — `T*` in the Casazza–Lynch letters the corpus uses; Christensen calls it `T` — a different object from the dual-frame reconstruction). Borrowing a vocabulary obligates you to its internal distinctions, not only its flavour.
- check: a collision across literatures renames only the colliding member: "projection" already names the idempotent `R∘M`, so the projector keeps its established name and the face takes its frame-theory one.
- check: let the inputs name the outputs (`Frame(basis, measure)` has the basis space as codomain and the measure space as domain), and keep the name that emerged as understanding deepened over the habitual one: it is the name that re-steers the next context. The secondary vocabulary (the FEM synonym beside the frame-theory primary) lives in the docstring, since consumers come from several literatures.
- check: generic machinery is named by STRUCTURE, never by a domain role: `Composite[Interior, Boundary]` over `System[…]`; the domain reading is the consumer's, carried by the specialization (`[R]` the user, 2026-07-10).
- **One word-order per family, the qualifier token always present and in the same place**, so that a grep over the shared token finds the whole set (`[R]` the user, 2026-06-22: "if we're perfectly and predictably consistent, you will know exactly what to search for with a grep, and so will anyone else"). check: before minting or renaming a member, grep the siblings, extract their word-order and token, and follow it exactly; an off-pattern sibling is renamed in the same change.
- tell: a member a grep over the family's token misses (`HarmonicMomentField` in the Domain-first `<Domain><Role>` family of `AngularFlux`, `ScalarFlux`, `BoundaryFlux`: renamed `HarmonicMomentFlux`, with `HarmonicMomentSourceSink` and `HarmonicMomentProjection` completing the matched set, and a role-first proposal refused for introducing a third word-order); a name patched with an adjective; a borrowed word whose literature-internal distinction the object does not honour.

## Values are loaded or computed, never typed

Pythonic code: dataclasses, type hints, scipy. Every number a module or a test carries is produced programmatically: a reference eigenvalue or coefficient comes from the derivation that defines it (`derivations/`) or the data file that carries it, never from a table by hand. check: a literal with five or more significant digits outside `derivations/` names its source or is replaced by a call. tell: a literal that matches a published table; a test pinning a number whose only provenance is a comment.

## A bare `assert` outside a collected test module is not a contract — the canonical runner strips it

`python -O -m pytest` is canonical; `-O` sets `__debug__ = False` and removes every `assert` at compile time. The scope: pytest rewrites the `assert`s of a COLLECTED test module, so those survive; a bare `assert` anywhere else, in production code, a test helper, a fixture, a `conftest` or a generator, is compiled out. A contract written as a bare `assert` there **does not run in the suite that matters**: production ships accepting the input the assert refuses, and a helper's gate passes whatever it is given.

- check: `grep -n "^\s*assert " orpheus/` and sort the hits. **Type-narrowing** (`assert x is not None` for pyright) may stay: it was never the guard. A **numerical / domain / admission contract** (a tolerance, an invariant, a shape law) **MUST be a real `raise`**, modelled on the nearest admission guard (`_assert_alpha_dome_closes`) so the vocabulary stays greppable.
- check, in tests: an assertion outside a collected module is `np.testing.assert_*` or a `raise`; a test that exercises a production `assert` on purpose (a Layer-3 check) carries the `skipif` pair of `vv-testing`.
- **Prove it, don't argue it.** Run the guard's own arithmetic on a deliberately-bad input under `python` and under `python -O`. tell: it returns instead of raising, so the contract is inert. [case](../evidence/coding-standards.md#2026-08-12-alpha-dome-assert)
- **Converting one is a retirement** (the `retirement-audit` skill); tests pin the **shortest distinctive fragment of the OLD assert's message**: grep that, not your new wording.
- Cardinal Rule 2 first, then the guard: the founding recursion had **three copies**, which is why the contract could live on one arm only.

## A guard is elegance debt — tag it, and name what retires it

A runtime guard (`require_member`, `admit_composite`, a typed refusal on an alien carrier) is a **signal that the architecture failed** to make the mistake unspellable (`coding-elegance` Pattern 4). A legitimate protection *today*, not the target state: **the ultimate state is to not need the guard** ([R] user ruling, 2026-09-07: [case](../evidence/coding-standards.md#2026-09-07-guard-is-debt)).

- **Every guard that lands carries a greppable marker in its docstring:** the token **`ELEGANCE-DEBT[guard]`**, the issue number, and ONE sentence naming the structural change that makes the guarded mistake unspellable (e.g. *"retires when B is bound on its own trace end"*). check: `grep -rn "ELEGANCE-DEBT" orpheus/` is the debt ledger.
- The issue is filed **with the carve that lands the guard**, never before (a guard without its retirement plan is an unpriced debt; a plan without its guard is a promise). The step landing the structural change deletes guard AND tag in the same commit, and the mutation battery must show the mistake is now unspellable, not merely refused.
- tell: the guard's docstring justifies the *check* rather than naming the *shape that would make the check unnecessary*.

**A guard is debt, or it is a declared scope boundary.** A guard stands where the machinery that would derive the refused property does not exist. When that machinery is in scope, the guard is debt, as above. When a user ruling has decided not to build it yet, the guard marks the edge of what the code derives, and it moves only when the machinery is built ([R] user ruling, 2026-09-22; the model is the orbit catalogue, "deferred, not refused"). Either way the machinery is **named**, so that building it can be ruled on.

- **A boundary guard carries `SCOPE-BOUNDARY[guard]`** in its docstring, with `machinery:` (what would derive the refused property: "the Gröbner-basis orbit-space engine"), `ruling:` (who ruled, when, and where it is recorded) and `revisit:` (what would overturn it) within three lines. The same ledger test reads both tokens.
- **A boundary is well-formed** when it is single (one table and one door at one site), declared (its refusal says it is a scope edge and lists what lies inside), seeded (each entry is the missing machinery's output in the machinery's own data model, so building the machinery develops the entries instead of replacing them), and verified as if derived (each entry is tested against the law a derivation would have established).
- **An input boundary** is the other legitimate guard: data entering from outside (a data file, a user's deck) cannot be typed at its source, and is parsed into types once, at the edge, never re-checked inside.
- check: for any guard, name the machinery that would make it unnecessary. It exists or is in scope: debt. A ruling defers it: boundary. Neither: an undeclared boundary, a question for the user naming that machinery, never a verdict.
- tell: a boundary that decays into debt: an entry added outside the machinery's data model, a second catalogue or door, a consumer that routes around the door.

## Retire as you go — the audit is the `retirement-audit` skill

Superseded code is noise that invites extending the wrong path. **Retirement is a first-class deliverable.** Every refactor introducing a better pattern MUST retire its predecessor; shims live **one merge cycle only**; never keep backward-compat unless the user explicitly authorizes it. The audit is its own numbered substep, with a `file:line` retirement list.

- check: before calling any delete, rename or re-home done, run the numbered
  audit of the `retirement-audit` skill (qa, archivist, method-implementer,
  elegance-enforcer and test-architect hold it preloaded; the main agent loads
  it; a sub-agent's skills are its preload, since no Key agent lists the
  `Skill` tool): the three searches, the
  surfaces a symbol grep cannot reach, the migration of tests and markers, what
  the retirement does to the surviving gates. Its item B.4 has no other catcher:
  a dead docstring reference produces no build warning at any severity, and
  only `dead_references` reads that surface.
- tell: a retirement with no `file:line` list; a delete-only diff; a green
  build offered as the docs check.
