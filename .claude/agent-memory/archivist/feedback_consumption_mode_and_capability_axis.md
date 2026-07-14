---
name: consumption-mode-and-capability-axis
description: Two-page recipe for documenting a NEW consumption mode / capability axis added to an operator algebra (solve/apply/ASSEMBLE) — the per-method "third mode" narrative page + the numerics "new three-layer axis" page, tied by a bidirectional cross-ref.
metadata:
  type: feedback
---

When a refactor adds a THIRD way to consume one operator (beside
solve/apply — here **assemble**: emit the same per-cell coefficients as
`(row,col,value)`), the documentation splits across TWO pages with a
bidirectional `:ref:` seam. Instance: stencil-assembly 2b
(`refactor/spatial-promotion-assembly`, #272/#284/#282) — `loss_representations.rst`
§`loss-rep-three-modes` + `operator_algebra.rst` §`operator-algebra-assembly-axis`.

## Page 1 — the per-method "third consumption mode" narrative (the domain page)

Anchor on the page's EXISTING "canonical operations" framing. The native-frame
section already had a `.. list-table::` of *two* operations (SOLVE=sweep /
APPLY=matvec) as "the same operator viewed two ways" → **add the third row**
(ASSEMBLE) + flip "two" → "three", forward-ref the new section. Then the section
(H1) with these H2 subsections, in order:

1. **Intro + a 3-row mode table** (Mode / Consumer / Walk) — copy the carrier
   module docstring's table verbatim (it IS the algebra-of-record).
2. **The one-source guardrail** (the Phase-F twin-path lesson): assemble MUST
   read the coefficient source solve/apply read — never a second stencil
   spelling. Consequence: a sign flip in the shared kernel moves all three modes
   together (teeth bite BY CONSTRUCTION). "A leaf whose apply is monkeypatched to
   RAISE still assembles" is the crisp proof they share coefficients not call paths.
3. **Coefficient extraction by kernel probing** — the load-bearing WHY the teeth
   hold: the residual/trace maps are AFFINE (`is_linear`), so a UNIT probe reads a
   coefficient block EXACTLY (an algebraic identity, not a fit). Write the two
   affine relations (`r = Āψ̄ − ΣC_bψ_in − Q̄`; `ψ_out_a = T_aψ̄ + Σα_ab ψ_in_b`)
   + the probe recipe (cell probe `e_m⊗1` → column m of Ā/T_a; face probe `e_f⊗1`
   → column f of −C_b/α_ab). Show how the SAME extraction serves every closure
   (DD α=−1 dense chains vs LD α=0 block-tridiagonal).
4. **The walk-order triangularity theorem** — the object-level "sweep IS forward
   substitution": the symbolic walk carries each in-flight face's ROW BLOCK; each
   cell's rows depend only on already-emitted upstream rows → (block-)lower-
   triangular in walk order. THIS is the #284 discharge as a matrix fact.
5. **Frame conjugation** (multi-moment closures) — the sweep-frame → global-frame
   involution as a similarity `M_global = Φ M_sweep Φ`; tie it to the ERR-061 root
   cause (lift the slope moment before a frame-agnostic consumer sees it).
6. **Zero-inflow posing** — the bulk block is emitted at zero trace; the trace
   coupling A_bs is a SEPARATE composite block; every consumer poses at zero inflow.
7. **Scope + the obstruction** — Cartesian-only, with the curvilinear obstruction
   named as a walk-order BACK EDGE (#282), and the characterization gate framed
   POSITIVELY (assert the defect present, never xfail — the doc-level tripwire that
   reddens when the fix lands).
8. **Verification** — the gate families (G1 object-level matvec equiv @ ~6e-16,
   never a scalar functional (vv Mode 12); G2 LAPACK solve_triangular ≡ sweep =
   the #284 discharge, triangularity leg an EXACT structural zero; one-source teeth).

## Page 2 — the numerics "new capability axis" (the operator-algebra page)

Place right after the existing Mat-functor / materialising-functor section (the
assembly axis is its natural extension). Cover, as H2s:

- **Two Mat-functor realizations** — a `.. list-table::` contrasting apply-to-basis
  PROBING (total, dense, O(n) applications; the RETAINED `_as_matrix_by_probing`)
  vs structural EMISSION (partial — only where a stencil exists — O(nnz)). Frame the
  carrier as a scipy.sparse SERIALIZATION, NOT a new COO-builder algebra (that would
  twin the operator algebra one layer down).
- **The three-layer surface** — minted EXACTLY like inverse/adjoint (predicate
  `is_X` / narrowing Protocol `SupportsX` (NOT runtime_checkable) / `TypeGuard`
  bridge `x()` / eager `MissingX(TypeError)` refusal). Cross-ref the existing
  three-layer section's anchor (the section may only have a DIFFERENTLY-NAMED
  anchor — grep `_.*:` above the header; here `capability-set-semantics`, NOT the
  slugified title — don't invent `the-three-layer-operator-surface`).
- **Composer homomorphism laws** — `[A+B]=[A]+[B]`, `[AB]=[A][B]`, `[αL]=α[L]`
  (CSR add / matmul / scalar); `is_X` = conjunction over legs; the unbuilt
  member (TensorProduct→kron) DEFERRED with no consumer (Cardinal Rule 2).
- **The delegation ruling** — `as_matrix` keeps dense semantics and DELEGATES to
  densified emission when `assemblable`, probing stays the fallback; no-op for
  every pre-landing call site (bit-safety by construction).
- **The anti-tautology discipline** — once as_matrix delegates, a probed≡assembled
  gate MUST force the RETAINED probing pathway (else it compares assembly with
  itself — vacuous). This IS the fuller-view-oracle rule; one permanent pin per
  family.
- **First production consumer + the Mode-11 sentinel** — the consumer that now
  runs the new path "automatically" is EXACTLY where a coverage gate goes vacuous
  (vv Mode 11 — green twin that never executes the rewired line); pin it with a
  SENTINEL (a counter monkeypatched onto the rewired reader) asserting the new
  path fired, not just that values match. (Here: the diffusion resolvent
  `MatrixInverseOperator(FlattenedOperator(A,template))`; the sentinel counts
  composite-dimension probes and asserts ZERO — allowing the tiny `ng` law-probes
  the boundary emitter legitimately makes.)

## Changelog + L-007 landed-seam

- Add a dated Development-history row to EVERY page that has the table AND lands
  content (here BOTH operator_algebra AND discrete_ordinates — even though the
  task named only two pages; Cardinal Rule 3 keeps the changelog complete where
  the section lands). Match each page's table format (`:widths:`, in-dev vs dated
  Where column).
- A page with only PROSE history (no dated table) gets NO row — instead reconcile
  its "deferred extension point" subsection: here loss_representations had
  envisioned this exact mode as a deferred "ExplicitMatrix" FIFTH SCHEDULE; it
  shipped as a third MODE (not a schedule). Classic L-007 landed-seam-different-
  shape: preserve the WHY, flip tense ("was envisioned"), tell "envisioned X /
  shipped Y", cross-link. The anticipated "direct-solve cross-check" use case WAS
  the G2 gate — say so.

## Traps hit

- Em-dash `—` in a title is 1 code point but the underline still came up 1 SHORT
  (67 vs 68) — size with `len(title)` in python (L-009), the `-W` build catches it.
- The repoint-target eq-labels were FORMER ORPHANS (documented, 0 tests): pointing
  the repointed tests at them KILLED the phantom AND covered two orphans in one
  move (see lessons.md L-003 phantom-repoint bullet).
