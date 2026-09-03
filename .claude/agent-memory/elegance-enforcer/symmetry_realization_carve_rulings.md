---
name: symmetry-realization-carve-rulings
description: #434 R1 (group-as-realization) elegance review — the tag→realization functor shape, where the invariant tier lost its guard, and the reusable probes.
metadata:
  type: project
---

# #434 R1 — "every question about a group is computed from its realization"

Reviewed 2026-09-02/03 on `fix/angular-phantom-support` (uncommitted carve in
`orpheus/numerics/symmetry.py`). Rulings that transfer to R4/R2/R3 and to any
tag→structure carve.

**The shape that worked, and is worth copying.** A value type keeps a **TAG as its
IDENTITY** (`__eq__`/`__hash__`, naming, constructors — ndarray-bearing structure
cannot hash) and derives its **STRUCTURE** through ONE memoised functor
(`_realize(tag) -> Realization`). `[M]` 68 → 26 `isinstance` sites; a 109-line
`_contains` → a 5-line `Realization.contains`. Measure the carve with an AST
`isinstance`-by-owner census before/after — it grades the claim in one command.

**⭐ The recurring defect this carve produced, and the one to grep for on every
tag→structure carve: THE GUARD DOES NOT MOVE WITH THE MATHEMATICS.** Every *tag*
class had a `__post_init__` (and the façade gained one); the two NEW value types
that now carry the whole theorem (`IdentityComponent`, `Realization`) shipped with
none — both public in `__all__`. `[M]` constructible: `IdentityComponent((zeros,))`
→ `dim 1`, `.axis = [nan nan nan]` (RuntimeWarning, no raise);
`IdentityComponent((X,X,X))` → `dim 3` ⟹ `contains_element` True for EVERY proper
motion; `Realization(IdentityComponent(()), ())` → a "group" whose
`contains_element(e)` is **False**. Root cause: **`dim` is spelled `len(generators)`
where the mathematics means RANK.** ⟹ when a carve moves the mathematics to a new
tier, ask *which `__post_init__` moved with it* — the answer is usually "none",
because the guards were written against the old tier and nothing prompts a re-look.

**A three-way `dim` branch on {0, 1, 3} is ESSENTIAL, not a tag in disguise** —
so(3) has no 2-dimensional subalgebra, so a new group cannot force a new arm (the
skill's own discriminator). PASS it. But it is sound only if `dim` is validated;
the two findings are the same finding read from both ends.

**Reusable probes (`scratch/_r1_elegance_probe*.py`).**
- 1: AST `isinstance`-by-owner census, new file vs the pristine copy.
- 5/15: illegal-state constructibility on the new value types + a candidate
  `__post_init__` validated against **all 31 shipped members (0 refusals)** AND the
  four illegal states (all refused). Always validate a proposed guard both ways.
- 9/14: two-spelling coextensiveness at a stated denominator (**961 of 961** ordered
  pairs) — the leg-2 downgrade from VIOLATION to CONCERN, measured not asserted.
- 12: `ast` walk for `Attribute.attr == "_tag"` — finds every private-field reach
  from module scope (here exactly 1, in `orbit_certificate`).

**Two verdict lessons.**
- A **façade→leaf rename that DROPS the load-bearing qualifier** (`generic_images`
  → `Realization.images`) is a real naming finding: "generic" *is* the vv#13
  content, and one grep should find both ends.
- A retirement can **widen a capability silently**: `Dinfh.generic_images` went
  `NotImplementedError` → 24 arrays because the refusal had been an artifact of the
  tag dispatch. The tests were re-keyed; the production docstring was not. ⟹ after a
  dispatch→computation carve, re-run every *refusal* the old dispatch spelled and
  ask whether the new body still refuses.

**Name collision to watch (reported as CONCERN):** `SubgroupOfO3.realization`
(component + coset reps) now sits beside the pre-existing `Quotient.realization`
(the manifold an orbit space IS, `[M]` ~15 production reads in `manifold.py` /
`measure.py` / `descent.py`). Two concepts, one word, one subsystem. R2 puts both in
one frame.
