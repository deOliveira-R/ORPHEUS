---
name: double-category-architecture-insight
description: Recipe for documenting a deep CONCEPTUAL architecture insight (a categorical/structural framing of an already-built type system) into an existing theory page — grow the seeded section into the structural treatment + the HONEST design rationale (what-was-explored-and-why-the-ideal-is-impossible), cross-ref not duplicate, all new eq-labels are vv-status:documented representational identities. Instance: the (Rep × Role) carrier grid as a DOUBLE CATEGORY in operator_algebra.rst (P4.5, #268/#261).
metadata:
  type: feedback
---

# Documenting a conceptual architecture insight (the categorical-framing-of-an-existing-type-system task)

A recurring task TYPE: the codebase already SHIPS a type system; the user
wants the DEEP structural/categorical framing of it archived (Cardinal
Rule 3 — "Sphinx is the LLM's brain"), with the emphasis on the HONEST
design rationale (the what-was-explored-and-why-the-ideal-can't-be-built
section). NOT an in-flight implementation status — frame as the
realized/intended design.

**Why:** the user explicitly asked for the "what we tried and why it
fails" content — the impossibility argument is the highest-archival-value
part (it stops a future session re-attempting the obvious-but-impossible
collapse and hitting the wall the hard way). This is the success-story
analogue of the close-out arc's structural-obstruction step.

**How to apply** (the 6-move recipe, instance: carrier grid as double
category, `operator_algebra.rst`, P4.5):

1. **Grow the seeded section, don't bolt on a new page.** The page
   already had a `scattering-carrier-grid` 2×2 sub-grid (seeded P4f). The
   insight LIFTS it to the full structure. Insert new H3 subsections
   bracketing the existing material: the structural lift right after the
   seed, the census + rationale as the capstone before the existing
   "Deferred relocation". Cross-ref the existing affine-torsor /
   frame-seam / heteromorphic-apply sections — do NOT re-derive them.

2. **The structural framing as a `.. list-table::` mapping
   categorical-part → in-the-code → consequence-for-the-code.** A double
   category: Objects=cells, Horizontal=Rep-morphisms (role-generic),
   Vertical=Role-morphisms (rep-generic), the 2-cell=the conjugated
   operator. The THIRD column ("consequence for the code") is what makes
   it a real architecture doc not abstract nonsense — every categorical
   part must cash out as a code fact (which generic, which 0-ULP gate).

3. **The interchange/coherence law IS a theorem the code already pins.**
   Find the existing bit-identical (0-ULP `np.array_equal`, NOT allclose)
   crosscheck and RE-FRAME it as the coherence witness. The 0-ULP (vs a
   tolerance) is the LOAD-BEARING point: a tolerance would admit two
   different reduction trees agreeing to round-off → the square wouldn't
   commute; 0 ULP says it commutes exactly = the mark of a real 2-cell.

4. **The census as a stub-columns `.. list-table::` (Rep rows × Role
   cols) + the PRINCIPLED HOLES.** An empty cell that is principled
   ("do not mint this") is as load-bearing as a populated one. Document
   each absence's REASON (the `(Moment,Residual)` hole: moment space is
   never the subject of a balance → no `from_balance` consumer). The
   role-axis "asymmetry" (some roles are mixins, some bare) is the
   type-vs-property rule WORKING (mint a role-object only where a
   non-identity morphism lives), NOT a defect — say so explicitly.

5. **The impossibility argument as a numbered-obstruction
   `.. list-table::` (# / obstruction / why-fatal).** Each obstruction
   fatal on its own. The conclusion: the flat normal form is the UNIQUE
   principled form, NOT a compromise. Close with a `.. note::` "this is a
   closed exploration, recorded so it is not re-opened" + pointer to the
   cross-domain-attacker memo that produced the verdict. This is the
   success-story version of the close-out's structural-obstruction step.

6. **Realization-status via a Development-history changelog line, NOT
   "P4.5 will…".** If the page lacks a Development-history section, ADD
   one (mirror the SN page format: `.. list-table::` When/Milestone/Issue/
   Where; `in dev (date) + branch` for unmerged). Frame the landed part
   as landed, the not-yet part honestly (here: two-param operator type
   LANDED on branch; the `@overload`-confession retirement NOT yet done,
   tracked #261). Trust-git caveat in the intro paragraph.

## The accuracy discipline that mattered most (L-001 in action)

The brief described P4.5 as if PENDING ("P4.5 will…"), but **grepping the
live branch showed `LinearOperator[Domain, Codomain]` was ALREADY in the
code** (operator.py:282/92-93, the double-category comment already in the
source), while the `@overload` confessions were STILL live (W-F not
landed). Reading the live code (not the plan's tense) is what let me frame
the realization status honestly. ALWAYS grep the live branch state before
writing "will/realized" — a surgical-carve plan's tense lags the code.

## V&V discipline (L-004)

All new eq-labels are STRUCTURAL/REPRESENTATIONAL identities (the grid-cell
layout, the interchange-coherence identity, the operator-typing identity)
→ each gets `.. vv-status: <label> documented` with a rationale comment
naming the verifiable content (the 0-ULP crosscheck, the foundation leaf
tests, the `assert_type` pins). They land in the matrix's "Documented-only"
bucket — verified by reading `docs/verification/matrix.rst` after the build.
NONE is a solver/eigenvalue/flux claim → zero prohibited V&V phrasing.

## Build/cross-ref facts confirmed this task

- `orpheus.numerics.operator` is NOT automodule'd → `:class:`LinearOperator``
  and `:data:`Domain`/`Codomain`` render PLAIN TEXT by page convention
  (the page already uses the `:class:` form 18× — matching it is correct,
  NOT half-surfacing). L-002.
- `:class:!HarmonicMomentResidual` (the `!` prefix) renders the
  deliberately-ABSENT class as a plain literal with no broken xref — the
  right spelling for a principled-hole reference.
- Baseline AND post-edit `-E` build: 0 WARNING / 0 ERROR / 0 CRITICAL,
  EXIT=0. The 10 "has no matching equation node" lines are INFO (verifies
  section-anchor registrations), not warnings.

Cross-refs: [[lessons]] (L-001 live-code, L-002 plain-text xrefs, L-004
vv-status), [[feedback-carrier-grid-typed-seam-layering]] (the P4 grid
seed this grew from), AGENT.md Close-Out arc (the structural-obstruction
step is the template for move 5).
