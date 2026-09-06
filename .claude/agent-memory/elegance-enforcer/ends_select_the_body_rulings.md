---
name: ends-select-the-body-rulings
description: CS4c step 5 (per-call carrier dispatch → one construction-selected body) — the four probes that find what a dispatch-to-admission carve loses, and the shapes it produces
metadata:
  type: project
---

# CS4c step 5 — retiring per-call carrier dispatch for construction-selected bodies

Reviewed 2026-09-05, branch `refactor/cs4c-step5-construction-selected-bodies`
(uncommitted production diff + `orpheus/transport/operators/{lift,angular_lift}.py`).
Verdict CONCERNS RAISED — 0 blocking, 11 should-fix, 8 nits. Design record
`.claude/plans/cs4c_binding_design.md` §19–§20 (rulings R-1..R-6).

**Why:** this carve shape (N `singledispatchmethod` tables + M `isinstance` carrier
arms → one body selected at construction from the operator's two ends) is the SN
campaign's recurring move, and it recurs with the SAME four losses. **How to apply:**
run the four probes below on any dispatch→admission carve before writing findings.

## ⭐ The four probes (run these; do not read for them)

**P1 — does the new admission see the ROLE?** A carve that replaces carrier-CLASS
dispatch with SPACE-CONTENT admission loses the role dimension, because in this
codebase **role is class identity and the space is SHARED across the pair** (the S2a
doctrine, `Field._check_partner`'s own docstring). `[M]` 2026-09-05:
`AngularFlux.zeros(s).space is AngularSourceSink.zeros(s).space` → **True**. So a
SOURCE fed where a FLUX is required passes the space check and dies later as a bare
`AttributeError` from a `cast()` that lied — where the retired dispatcher raised a
`TypeError` naming the operator. Probe: build the operand's role PARTNER on the same
space, feed it, read the exception TYPE. Fix shape: the construction-selected strategy
declares the carrier class it reads; the admission narrows once (not per call, so no
dispatch smell and no AST-gate trip).

**P2 — does a COMPOSITE admission admit both blocks?** `admit_composite` checked
`x.interior.space` and nothing about `x.boundary`. `[M]` fed a composite with a correct
bulk and a foreign mesh's trace: ADMITTED, and the emitted `FullField`'s flat dimension
(576+1920) violated the operator's own declared `codomain` (576+672). Probe: build the
operand from mesh A's bulk and mesh B's trace. Note the tree carries BOTH conventions —
`lift_bulk_action` rides the OPERAND's trace (CS4b S4), `RadialCharacteristicEmission`
rides the OPERATOR's.

**P3 — is the guard the code declares actually LIVE?** Run a positive control, do not
read the comment. `[M]` `AngularLift.__post_init__`'s per-END energy admission is INERT
on every shipped binding (a subclass returning `data_ng = 999` CONSTRUCTS) because
`FullFieldSpace.from_blocks` always yields `axes=None` and the guard reads
`space.axes`. Consequence: an abstract member of the subclass contract exists to feed a
check that cannot fire, while the REAL guard is an eagerly-bound sub-operator's own
admission. Probe: subclass, lie about the datum extent, construct.

**P4 — does the carve hand-roll a structure the tree already memoises?** `[M]`
`FunctionSpace.of_axes(*interior.axes[1:])` (drop axis 0 by POSITION) vs
`interior.retraction("angular").codomain` (drop by LABEL, frame-induced, memoised, and
the space `integrate_angular()` actually returns): **content-equal on 2 of 2 chart
families, different objects, 39× slower (12.97 µs vs 0.33 µs)** — and read per-apply on
the windowed hot path. Probe: `timeit` both and compare `==` / `is`.

## The two shapes worth recognising again

**⭐⭐ A string-keyed branch whose two arms are PROVABLY EQUAL is the worst kind, because
its miscompilation is numerically invisible.** `transfer.py:442`
`self._end.__name__ == "_MomentEnd"` — a module-private class named as a string literal
from another module (census: 6 references, 5 real + this 1 string; no import). `[M]` the
two arms it selects are **bit-identical** (`max|Δ| = 0.000e+00` on a random anisotropic
moment operand). So a rename flips the branch, every gate stays green, and what dies is
the *deliberately-kept surface* the other arm is the sole consumer of (here R-5's typed
`HarmonicMomentFlux` arm + the minted `source_reconstruction` face) — which the next
dead-code audit then retires as unconsumed. ⟹ when two routes are ruled equal-and-both-kept,
the SELECTOR between them needs a construction-time `KeyError`, not a string, and the
mutation battery needs an arm that mutates the selector (it will read GREEN today —
that green IS the finding).

**A hoist to a new base is when a PRE-EXISTING half-admission must be completed.** The
face-binding check (`flux_analysis.domain == emitted`, `source_reconstruction.codomain
== emitted`) was copied verbatim from `TransferOperator` into `AngularLift` — including
its gap: the two faces' MOMENT ends are never compared, so `[M]` a mixed-frame pair
(M at L=1, R at L=0, moment ends `(2,3,…)` vs `(1,1,…)`) CONSTRUCTS. The defect that
gap admits is **the one this family already shipped once** (#426 step 2: "two mint bodies
of one recipe, one of which minted at L=0"). Grade it should-fix (pre-existing ⟹ not
blocking by the diff-boundary rule) and cite `coding-standards`' clean-before-extend.

## Recurring local smells this carve produced (short list)

- **Two spellings of "derived once at construction" in one carve**: a declared
  `field(init=False, repr=False)` assigned in `__post_init__` (good) vs a
  `cached_property` force-evaluated by a **bare expression statement** (`self.x` on its
  own line — ruff B018, and the line a cleanup pass deletes; here it carried the only
  live ng guard).
- **Three spellings of one predicate in one file**: `isinstance(d, FullFieldSpace)` vs
  `isinstance(d, FullFieldSpace) and d.interior_space is not None`, twice. `[M]` they
  genuinely disagree (apply takes the composite arm, assemble the plain arm, on the same
  object) — both loud, so nit-grade, but invisible in review.
- **The same two properties byte-identical on two sibling classes** (`_domain_interior` /
  `_codomain_interior` on both lift classes). ⚠ Do NOT answer this by putting the parse on
  the universal `BoundOperator` — the plain-bound leaves would inherit a method that
  refuses on every instance (harmful partial method). A `CompositeBound` mixin is the
  home.
- **A `Literal[...]` two-valued tag beside a real `Enum` two-valued tag in ONE new module.**
- **A private static reached across a tier by 4 of its 5 call sites** — the underscore
  misdescribes the audience.
- **A docstring universal with no denominator**: *"the family's ONLY two carrier parses"*
  while the package still carries three more by ruling. Same defect class as
  `plan-authoring` §2's quantifier clause, wearing a module docstring.

## What this carve got RIGHT (reinforce; verified, not read)

- `assemble()` reproduces `apply()` **bit-identically**; `as_matrix` == `assemble`
  bit-identically; `BulkLift.assemble()` is index-identity on the bulk block with **all
  trace rows and columns zero**. An assembly-embedding verb is cheap to verify — always do.
- The `RolePair.__init_subclass__` bijection (one declaration registers both directions,
  a second source/sink for one flux is refused at import) is Pattern 4 at the right tier,
  and it is what let the per-family carrier parse die with no replacement branch.
- `on_moment_domain()` through `dataclasses.replace` — every admission re-runs and the
  construction selection re-lands; idempotent on an already-moment binding.
- A STRICTER array admission (exact shape) is a win, but verify the two consumers a
  strictness flip could redden — here the `Fold ∘ K ∘ integrate` ray row and the assembly
  impulse probes, both `(ng, nx)` exactly.
