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

---

# #434 R4 — "the lift is a derivation output, and an orbit space's dimension is a theorem"

Reviewed 2026-09-03, same branch, uncommitted carve in `manifold.py` (bulk) +
`symmetry.py` + `basis/descent.py`. Probes `scratch/_r4_elegance_probe{1..8}.py`.

**⭐⭐ The flagship, and it generalises past this campaign: a DERIVED field excluded
from `__eq__` + a `functools.cache` keyed on the OWNING object = call-order-dependent
wrong answers.** `Quotient.lift_codomain` shipped `field(compare=False)`, so
`replace(entry, lift_codomain=SPHERE) == entry` AND `hash` equal; `barycentre` is
`@functools.cache`d on the `Quotient`, so `[M]` `barycentre(entry).codomain` reads
`S^2` or `D^3` **depending on which of the two entries was asked first** — the
honest entry gets poisoned by the forged one. The field's own docstring says its
purpose is "the codomain a consumer reads to learn which it was handed". ⟹ **when a
diff adds a field whose stated job is to let a consumer DISCRIMINATE, check (a) that
`__eq__` can tell the two apart and (b) every memo keyed on the owner.** The stated
exclusion reason ("derived from `(base, by)`") proves too much: `realization` is
derived the same way and IS compared; the real reason the sibling callables are
excluded is that *a function has no value equality*, which a `Manifold` has.

**⭐ The R1 lesson recurred one carve later, in the same `__post_init__`.** R4 gated
`realization` (the dimension law) and left the two NEW required fields
(`lift_coordinates`, `lift_codomain`) with no clause — no width check against the
base's ambient space, though `_act_through` hands the lift's output to
`RigidMotion.on_points`. `[M]` `replace(entry, lift_codomain=Ball(2))` constructs
while the lift still emits width 3.

**⭐ A guard resting on an unstated "generic point" is an assumption where the
mathematics offers a computation.** `_generic_orbit_dimension` takes the rank at
`_GENERIC_SPHERE_PROBE[0]`. `[M]` poison row 0 with `e_z` and the dimension law
**inverts**: the true `S^2/O(2)_z` is REFUSED and the `Ball(2)` forgery it exists to
refuse is ADMITTED. The probe constant is a *seeded RNG* whose docstring names only
`descending_slots` as its consumer, and it has already been re-generated once in this
codebase's history. `dim(generic orbit) = max_p dim(H·p)` (orbit dimension is upper
semicontinuous) ⟹ **`max` over the probe set is both the theorem and the fix.**
Coextensive today (`[M]` 9 of 9 probe rows agree on all 7 continuous groups) ⟹ graded
CONCERN, not VIOLATION, under the three-leg standard.

**Native place, with the strongest form of the argument: the new function is the
GENERAL form of a method the target type already has.** `manifold.py` open-codes
`rank[X p : X ∈ 𝔤]` by reaching `group.realization.component.generators` — three
levels into another module's value types, and it is exactly `IdentityComponent.fixes`
generalised (`fixes(pts) ⟺ orbit_dimension(pts) == 0`). It also re-asks, with the
numpy DEFAULT tolerance, the question `IdentityComponent.__post_init__` already asks
with `tol=_ELEMENT_ATOL`. ⟹ when a carve adds a computation about type T outside T,
grep T for the special case — if one exists as a method, the general form belongs
beside it and inherits its band for free.

**Withdrawn attacks worth not re-running.** `is_trivial` as a tag test — `[M]`
`SubgroupOfO3.__post_init__` normalises `Cn(1)` to the `Trivial` tag, so the tag IS
the identity (R1's ruling) and the test is as strong as a structural one. A
`CoordinateChart(select, embed)` value type — fails the concept-count test (one field
+ one class replaces three fields; `quotient_map`/`lift` gain a hop) and its benefit
is absorbed by one `__post_init__` clause. `_embedded_nodes` now projecting a fold's
representatives — `π∘g∘P_H = π∘g` for `g ∈ N(H)`, so the chart answers cannot move
(`[M]` 0 of 9925) and the position test becomes strictly more correct.
