---
name: cross-solver-unified-law-doc-architecture
description: Recipe for documenting ONE cross-cutting numerical discipline (a law) realized differently across N solver families — canonical derivation on the architectural-prototype page + short sibling spellings that cross-ref it, per-method convention taxonomy, the dead-by-design consistency-theorem close-out for a retired injection seam, and the "protocol contract a mode WOULD implement" honest-drift fix for a non-conforming family.
metadata:
  type: feedback
---

When a cross-cutting DISCIPLINE lands (one law, realized differently
across the solver families), do NOT duplicate the derivation on every
method page. The reusable doc-architecture:

**How to apply:** any time a single correctness LAW spans ≥2 solver
families (k-estimator unification #259/#291; a future DSA-consistency
law; the #259 substrate-unification follow-up). Instance:
`refactor/k-estimator-unification` — the unified k-estimator law "the
reported k IS the eigenvalue of the fixed-source map every method
scales only fission by 1/k through".

1. **Canonical derivation on the ARCHITECTURAL-PROTOTYPE page** (SN —
   it carries the most development history because the algebra was
   built there first). Full treatment: the posed problem, the balance
   identity by **divergence-telescoping** (sum the discrete cell
   balance over all cells; interior faces telescope by current
   continuity → only boundary faces survive = leakage L — mirror the
   diffusion page's `1ᵀ(C−S)=Σa` discipline, cross-ref its anchor), the
   leakage functional (net_current contraction + face-measure-matches-
   volume_measure table + 3-D fail-loud), the structural-zero design
   (WHY: bitwise preservation of every reflective anchor — terms vanish
   *structurally* not numerically on the cases that must not move), the
   scale bridge + trace-of-last-inner-solve contract, the convention
   fork with the k_old = k* + e(1−k*)/a algebra, a characterization
   table, gotchas, verification. New section anchor (`.. _X-estimator:`)
   + eq-label (`:label: X-keff-update`).

2. **Sibling method pages carry a SHORT spelling + rationale +
   cross-ref to canonical — NEVER a duplicate derivation.** Classify
   each family by its role in the unified law (a per-method convention
   TAXONOMY): **symptom** (SN — had the #291 leakage omission), **root
   flip** (MoC R7 — emission moved numerator→removal-denominator),
   **already-consistent** (CP — only the SUBSTRATE differs, unification
   is a follow-up), **degenerate** (homogeneous — L=0 direct-eig, law
   auto-satisfied). Each sibling keeps its verifies-target eq-label
   (L-003: rewrite body, preserve label name), adds a `.. note::` with
   its own convention fork + the principled-re-baseline value (MoC
   1.125→1.25) + a `:ref:` to the canonical section.

3. **Retired injection seam → "dead by DESIGN, not dead by being
   unwired," backed by a CONSISTENCY THEOREM.** When a kwargs/callable
   injection seam is retired to hardwired methods: (a) the theorem —
   apply `1ᵀ` (the `Σ`) to the converged fixed point `(A−S)ψ*=Fψ*/k*`
   ⟹ the hardwired ratio returns EXACTLY k*; every estimator consistent
   with the posed problem collapses to one value ⟹ a *different*
   injected estimator is by construction inconsistent ⟹ illegal states
   unrepresentable (coding-elegance Pattern 4); (b) the honest-`A.apply`
   contract — the theorem needs the TRUE net-removal action, so a
   stubbed adapter now yields a visible non-eigenvalue (the correct
   failure, surfaced by design). This is the success-story twin of the
   close-out's structural-obstruction step.

4. **Non-conforming family (MC) honest drift fix:** state the protocol
   member as "the contract a mode WOULD implement," NOT an existing
   member — FIRST grep `orpheus/<pkg>` to confirm it has no such member
   (MC had no `compute_keff`), then connect its NATIVE estimator (MC
   cycle weight-ratio) to the same law IN EXPECTATION, and spell the
   unified law it would have to satisfy. (The brief may mis-say the doc
   "describes an existing member" when the current text already said
   "cannot" — L-001: verify the live doc state, sharpen rather than
   invert.)

5. **Standard gates** (all fire here): keep verifies-target labels +
   rewrite body (L-003); `.. vv-status: X-keff-update documented` on the
   NEW definitional law-label (L-004 — a governing/definitional identity
   is NOT a solver-correctness claim; it files under "Documented-only"
   in the auto-regenerated matrix, rationale-comment names the
   consistency gate + notes the `verifies(...)` marker-wiring is a
   test-side follow-up when tests are off-limits); reproduce the
   LOAD-BEARING numbers analytically (L-001 — the k_old/k* fork algebra
   and the L/A overshoot ratio, not just transcribe); L-002 grep gate
   (0 `:func:`/`:meth:` xrefs to retired symbols; remaining mentions are
   intentional ``code`` literals) + HTML-resolution check for every
   cross-doc `:ref:` (renders plain-text when dangling, no `-W` warning);
   generated matrix auto-regens (L-008, verifies-target counts preserved).
