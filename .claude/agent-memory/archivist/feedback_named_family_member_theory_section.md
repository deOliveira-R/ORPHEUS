---
name: named-family-member-theory-section
description: Recipe for a NEW theory section documenting a NAMED member of an invariant-keyed operator family (the "name is a promise backed by a test" bar) — the GreenOperator/#226-step-4 instance; the name-earned invariant table, the iterative-member V&V framing, the ordering-ruling edge table, the promise-driven §18.A "what failed" content, and the OPEN→RESOLVED issue tense-correction.
metadata:
  type: feedback
---

# Documenting a NAMED member of an invariant-keyed operator family

The reusable doc SHAPE for a section that documents a new NAMED subclass in
a family where **a name is a promise backed by a test** (taxonomy §13 bar).
Instances (both `refactor/inverse-as-operator`): `GreenOperator` — the
first *iterative* inverse (#226 step 4, `operator_algebra.rst`
§`green-operator`); and `MatrixInverseOperator` — the dense DIRECT inverse
(#226 step 5, §`matrix-inverse-operator`), which ALSO promoted the universal
`as_matrix` functor and closed #285 for products. The 4 STEP-5 deltas below
generalise the recipe beyond a pure member section; the core 6 moves held.
Reusable for any future invariant-keyed family member.

**Rule.** When a carve mints a NAMED family member whose name is EARNED by a
distinguishing invariant, the theory section is NOT "here is class X". It is
a 6-move argument that the name is a promise the code keeps.

**Why:** the family is keyed by the mathematical OBJECT, not the algorithm
(a Richardson-Green and a Krylov-Green are the SAME Green — so
`KrylovInverseOperator` was rejected: "Krylov" names the realization axis).
The doc must make the object/invariant identity legible, or a fresh session
re-mints a per-algorithm sibling.

**How to apply — the 6 moves (each an h2 subsection):**

1. **The mathematical identity, LABELED + vv-status documented.** The
   core identity (Green: `(A−B)⁻¹ = Σ_k (A⁻¹B)^k A⁻¹`, the
   A-preconditioned splitting = multiple-scattering Neumann series) gets a
   `.. math:: :label:` + `.. vv-status: <label> documented` (L-004 —
   literature-transcribed definition, NOT a solver claim). Grep-collision-
   check the label first (L-003). Cite the pillar refs as the file's
   convention dictates — here `[LewisMiller1984]_`/`[AdamsLarsen2002]_` are
   defined elsewhere and resolve globally, so CITE-not-redefine (L-006);
   the adjacent section used plain-text, both are fine.
2. **"The name is earned" `.. list-table::`** of the distinguishing
   invariants, each row = *what it asserts* + *why the generic parent
   FAILS it* (G-Neumann: the splitting a bare `A⁻¹` has none; G-reciprocity;
   G-kernel folded into the δ_j anchor). This is the load-bearing
   name-justification. The success-story analogue of the close-out's
   structural-obstruction table.
3. **"Wraps the driver, re-implements nothing" (§11.2, Pattern 5)** — the
   member is a thin wrapper deriving its config from STRUCTURE (left-spine
   head → preconditioner via its own structure-keyed `.inverse()`; rest →
   negated gains), handing to the SAME shared driver the solvers consume.
   The flatten walks EXACT-type nodes only (subclass = structural leaf, so
   the MRO shadow survives) — verify against the live walk
   (`while type(node) is OperatorSum`), not memory.
4. **The ordering ruling `.. list-table::`** (the #261 canonical-order rule
   extended) — the N edges of "operand spelling selects the algorithm":
   fuses→direct (MRO shadow) / canonical→converges / legal-but-divergent→
   constructs-then-raises-LOUDLY / non-invertible-first→refused-at-
   construction. State the keying predicate's HONEST scope
   (`is_invertible = a.is_invertible` = "leading-term-preconditionable at
   this spelling", spec §18.B — NOT spelling-independent invertibility;
   lockstep `CAP_SOLVE` keeps the faithfulness keystone).
5. **The promise / "what was tried and failed" (Cardinal Rule 3 gotchas).**
   The §18.A delta IS the Rule-3 failed-approaches content: the driver's
   increment under-delivers by ρ/(1−ρ) (Signature 9), so the promise reads
   the TRUE residual and DRIVES a refinement loop (check-only would
   false-raise for every ρ>1/2). Include the adversarial FP shapes the
   divergence tooth found (increment=nan; denominator-overflow false
   convergence `res=finite/inf=0.0` onto a ~1e154 iterate) — this is
   exactly the "gotchas that aren't obvious from the code" the standard
   demands.
6. **Verification `.. list-table::`** — the L0 + L1 gates, the claim layer,
   the mutation count (12 bite + 2 designed-green controls). See the
   ITERATIVE-MEMBER V&V FRAMING below.

**Iterative-member V&V framing (lessons L-010, reinforced):** the FIRST
iterative member of a previously-all-EXACT family has NO bit-id twin →
correctness rests on structural-independence anchors ONLY (dense-LU +
Neumann expansion), and carries NEITHER an eigenvalue NOR an MMS claim
(an iterative sum inverse is neither an eigenvalue solver nor source-driven
— both pillars INAPPLICABLE). When the parent's `solve` is DEFINED as
`inverse().apply`, that equivalence is a TAUTOLOGY for this member — exclude
it as evidence EXPLICITLY. The name-earning invariant (G-Neumann), not
round-trip, is the correctness evidence.

**Euclidean transpose vs metric adjoint (L-010).** A "reciprocity"
invariant almost always uses the EUCLIDEAN transpose `Gᵀ` (built manually
over transposed operands, pinned by `⟨φ₂,Gφ₁⟩=⟨Gᵀφ₂,φ₁⟩`, no angular Gram),
NOT the metric Hilbert adjoint `.H = 𝒢⁻¹Gᵀ𝒢` (the #280 family). A taxonomy
table may aspirationally write `G.H`; what LANDED is the transpose. Write
the LANDED object; put the `.H`/#280 distinction in an `.. important::`.

**The OPEN→RESOLVED issue tense-correction (an L-007 instance on an
issue-DECISION, not a code-seam).** A carve often RESOLVES an issue a prior
doc passage framed as OPEN/deferred ("#285: structural vs per-leaf is the
open decision, lands when GreenOperator does"). Rewrite tense-correctly:
PAST-tense the interim state ("when step 3 shipped, `_seeded_inverse`
narrowed per-leaf…"), PRESENT the resolution ("step 4 resolved #285
STRUCTURAL — the mixin's abstract signature; pyright + ABCMeta enforce it"),
and NAME the residue for the next step ("composed `OperatorProduct.inverse()`
joins at step 5"). VERIFY the shape that shipped against the canonical
docstrings (here `SupportsSeededApply`/`_seeded_inverse`/`InverseWrapMixin`
are the canonical wording — read them, don't infer). The GitHub issue may
stay OPEN (residual scope) even though the core decision resolved — say
"resolved STRUCTURAL … #285's remaining scope" to capture both. Retitle the
section if its title carried "open" ("…open structural decision" → "…
structural resolution"; size the underline in code points, L-009). Also
sweep for convention stragglers a sibling missed (here an `L`→`A`
forward-operator name the L→A sweep skipped in the KrylovAcceleration
bullet — verify the target vs the live `def`, L-011).

**Homes:** the member section goes adjacent to its family's driver/factory
section (Green next to `inverse-application-driver`). The Development-history
CHANGELOG entry goes in the campaign's canonical changelog — here
`discrete_ordinates.rst`'s "Development history" (the SN page owns the #226
taxonomy changelog even though the member section lives in
`operator_algebra.rst`); reverse-chronological, above the prior step,
cross-doc `:ref:`<member-anchor>` (:doc:`operator_algebra`) — verify it
renders as an `<a href=…>` link, not plain-text (L-002 rendered-HTML audit).

**Verified against live code every claim** (L-001): the flatten walk, the
`is_invertible`/CAP_SOLVE lockstep, the refinement-loop apply body, the FP
shapes, the `_SolveBackedLeaf→_InvertibleForward` rename (grep-gate: gone
from `orpheus/`). Build `-E -W` EXIT 0, WARNING/ERROR/CRITICAL set unchanged
from baseline; matrix.rst auto-regenerated the 3 documented labels into the
vv-status bucket (leave as built, L-008).

**STEP-5 APPLIED (2026-07-02, `matrix-inverse-operator` §) — confirmations
+ 4 deltas.** The recipe held; the member section landed clean (`-E -W` EXIT
0, ZERO warnings — cleaner than the 1-warning baseline; 3 labeled eqs
`matrix-functor-out`/`matrix-inverse-materialise`/`matrix-inverse-direct-residual`
all `.. vv-status: … documented` → auto-regen'd into the "Documented-only"
matrix bucket). Deltas worth carrying forward:

1. **When the member carve ALSO promotes a universal BASE method, the section
   grows from 6 moves to ~8 subsections and documents BOTH.** Step 5 promoted
   `as_matrix` (the functor OUT of the operator category, `Op→Mat`, NOT an
   endofunctor — the taxonomy §2 "fourth arrow") alongside the leaf
   `MatrixInverseOperator`. The functor gets its OWN front subsections BEFORE
   the member ones: (a) functor-out framing + labeled apply-to-basis eq; (b)
   the resolution rule + TWO error classes (ill-posed `ValueError` ≠
   resource-refused `MatrixTooLarge` — class-DISCRIMINATED, never a loose
   `except (ValueError, RuntimeError)`); (c) "resource effect on a TOTAL
   functor → RuntimeError not ValueError, and NO `is_materializable` predicate"
   (§17 A2 — every operator HAS a matrix, so un-materializability is a resource
   fact, not a structural restriction like `is_invertible`). This
   materialize/serialize framing is reusable for any future `Op→Mat` method.
2. **The "supersedes the prior parenthetical" doc move (success-story analogue
   of a tombstone).** A universal-method promotion can make an OLDER
   name-earning claim FALSE — here §13's M-row "matrix-free CANNOT satisfy
   M-materialise (no `as_matrix`)" went false the moment `as_matrix` became
   universal (an iterative Green now also satisfies `[G][A]≈I`). The doc must
   EXPLICITLY say "**this supersedes X**" and RE-EARN the name on a SHARPER
   axis — the precision GRAIN (machine·cond vs driver-tol), with a CONTRAST
   `.. list-table::` (same `A` wrapped as the iterative sibling meets only the
   weaker grain). Not a retraction (no `.. note:: Retraction`) — a SHARPENING;
   the claim gets stronger, not withdrawn.
3. **Recipe move 4 (ordering-ruling edge table) is CONDITIONAL — it applies
   only if the member is a `.inverse()` FACTORY-DISPATCH target.** A
   DIRECT-CONSTRUCTION-ONLY member (no auto `.inverse()` routes to it — spec
   §35.3) has no "spelling-routes-the-algorithm" table; the VALUES-vs-STRUCTURE
   WITNESS table replaces it (the ctor reads VALUES not STRUCTURE → it inverts
   the leading-non-invertible sum the STRUCTURAL sibling REFUSES; the
   `(−S_ao)+D` ndarray analog of the FullField `(−S)+(L+C)`, since typed
   carriers are out of `as_matrix`'s ndarray scope — an honest-scope `.. note::`).
4. **A retire-private-helper + re-role-a-method step has a cross-page STALE
   CONSUMPTION-CLAIM blast radius (grep the WHOLE `docs/` tree, L-002/L-007).**
   Step 5 retired the private `_as_dense` AND re-documented `dense_per_material`
   (assembly-path → test-only storage ORACLE). BOTH `homogeneous.rst` (theory)
   AND `api/homogeneous.rst` had CLAIMED the K_iso assembly ran via
   `dense_per_material` — a Cardinal-Rule-1-wrong consumption claim (it never
   did; now definitively an oracle w/ ZERO production consumers). Fix: repoint
   the mechanical attribution to the REAL path (`as_matrix`), reframe
   `dense_per_material` as the oracle, but KEEP the `fission-source`
   verifies-target label + its math (L-003 — fix the attribution, not the
   claim). SAME-CUT sibling staleness rides along: the F-dyad `numpy.outer`
   claim (both homogeneous call sites flipped to `as_matrix` TOGETHER) —
   in-radius, fix it; the `numpy.linalg.solve`+`eig` eigenpair bullet
   (`direct_eigenvalue` internals, unchanged) — accurate, leave it.

Related: [[stub-to-rich-narrative-expansion]] (source-reading order),
[[operator-reification-retype-doc-pattern]] (the prior #226 step's carve
doc), and lessons L-004/L-007/L-010.
