---
name: orientation-axis-two-frames-doc
description: Recipe for documenting a 2×2-face operator unification where ORIENTATION (fwd↔adjoint) is the coherence axis and execution strategy is a NON-FREE third axis; plus the swap-law object-identity, the Euclidean-transpose landmine, and the verify-a-doc-fix-rider secondary-staleness trap.
metadata:
  type: feedback
---

Recipe for documenting an operator whose actions form a **2×2 of faces**
``{forward, transpose} × {solve, apply}`` that unify into **two frames,
each shared across ORIENTATION**, with execution strategy as a *non-free
third axis*. Success-resolution architecture doc (NOT the 9-step CLOSED
arc). Instance: #280 Phase 2.5e — the SN loss-operator walk unification.

**Why:** this shape recurs whenever a refactor adds the adjoint direction
to an already-forward-unified walk. The archival value is teaching the
next session that the coherence axis is ORIENTATION, not kernel — so the
forward solve-scan and apply-loop are NOT twins to merge.

**How to apply — the 8 moves** (loss_representations.rst
`loss-rep-orientation-two-frames`):

1. **The 2×2 faces `.. list-table::`** (solve/apply × forward/transpose).
   NAME the previously-empty cell (`sweep_transpose = (L+C)⁻ᵀ`).
2. **"Why it fragmented — the honest story."** The forward unification
   was REAL but scoped `{forward}×{multi-D}`; the adjoint + 1-D matvec
   arrived later as verbatim relocations, each standalone. Flip tenses,
   preserve the reasoning (L21 "achieved for multi-D forward, only
   CLAIMED for the 1-D walk + adjoint").
3. **The two frames, each shared across orientation.** apply-loop =
   `{loss_action, loss_action_transpose}`; solve-scan = `{sweep,
   sweep_transpose}`. List the fork points that orientation carries:
   reversed cell order + boundary **in↔out swap** + curvilinear Carlson
   **mirror routing** + the **second triangular factor** (`angular_adjoint`).
   Orientation is an OBJECT, pinned by an AST tripwire banning
   `is_adjoint/is_forward/is_transpose/is_reverse` (the one-walk
   `is_solve` tripwire's sibling).
4. **The orientation × kernel × execution taxonomy `.. list-table::`.**
   Coherence axis = ORIENTATION (free); execution `{scan, cell-loop,
   wavefront-graph}` is the NON-free third axis pinned by `(kernel,dim)`:
   solve→scan (inverting a recurrence is what scans are FOR — a per-cell
   loop loses vectorization), apply→cell-loop, multi-D→wavefront-graph.
5. **The Euclidean-transpose `.. admonition:: :class: warning`** (the #1
   landmine): the transpose is reversed-DAG + face in↔out swap over the
   SAME per-ordinate DAG — NOT μ-reversal, NOT the continuous adjoint.
   The physical Hilbert `.H` (`G⁻¹AᵀG`) rides ON TOP via the metrics; the
   sweep carries NO metric code.
6. **The deferral ledger** — every typed/loud `NotImplementedError` the
   unification must PRESERVE (multi-D Cartesian adjoint, LD-slab adjoint,
   ScheduledInvertibleOperator transpose, 2-D LD). Honest deferrals, not
   gaps.
7. **The swap law as OBJECT IDENTITY**, not a numerical coincidence:
   `A.H.inverse() ≡ A.inverse().H` because `_AdjointOperator.inverse()`
   returns `inner.inverse().H`. The metric adjoint-solve
   `G⁺·solve_transpose(G·b)` then falls out of the EXISTING
   `_AdjointOperator.apply` FOR FREE — no `_AdjointOperator.solve`, no
   metric in the sweep (a whole planned deliverable "dissolves"). Give
   the two-factor honest predicates (`is_adjointable`, `is_invertible`).
8. **The initial_guess / warm-start architecture split.** Direct exact
   inverses ACCEPT-AND-DROP `initial_guess`; the genuine warm start (the
   "3 vs 30 iterations" lever) lives at the ITERATION layer
   (`SourceIteration.solve` / `GreenOperator`). Made possible only
   because route (a) computed the seed DIRECTLY from source.

**Traps (verified 2026-07-05):**

- **ENDOMORPHIC metric = generic `G`, not `G_V`/`G_W`.** For a
  domain=codomain composite, use `G⁺·solve_transpose(G·b)` matching the
  code docstrings. Introducing `G_V`/`G_W` subscripts invites a
  subscript-ORDER error (the a3 plan's `G_W⁻¹·solveᵀ(G_V·b)` swaps them
  vs `_AdjointOperator.apply`'s `G_V⁻¹·applyᵀ(G_W·)`). Since V=W, drop
  them.
- **Page cross-ref convention is FILE-LOCAL and must be probed, not
  assumed.** `loss_representations.rst` renders ALL `:meth:`/`:class:`
  code-refs as plain `<code>` literals (0 py-links in the built HTML,
  INCLUDING pre-existing refs); `discrete_ordinates.rst` IS a linking
  page (62 py-links) yet its dev-history-TABLE code-refs still render as
  unresolved xrefs (same as the neighbouring rows). So: match the
  surrounding row/page style, write CORRECT dotted paths (verified vs
  LIVE code — they render plain either way, NO `-W` warning), and put
  the load-bearing navigation in `:ref:`/`:eq:` (those DO resolve
  cross-doc — grep-gate the built HTML for `href=".*#label"`).
- **V&V honesty:** a seed/closure re-pose re-baseline is justified
  STRUCTURALLY (the honest single-pass direct inverse — cite the cold
  residual collapse), NOT by angular accuracy. Pair the honest-scope
  note with a prophylactic `.. warning::`: the eigenvalue re-pose is
  certified by an angular-order N-sweep at fixed mesh, NEVER by MMS
  (Mode 7 — MMS is blind to the seed by construction). NEVER write
  "MMS verifies the eigenvalue".

**The "verify a doc-fix rider" secondary-staleness trap.** A rider task
names ONE flagged line (here: a `psi_half_seed : PsiHalfAngleSeed`
param doc naming a retired `CarlsonInwardSweep`). The implementation
phase usually already cleaned THAT exact line — but the stale
CONVENTION it encoded ("adopt φ½=0 for the matvec … converges under
fixed-point iteration") has a SECOND home in a sibling docstring (the
module docstring), findable ONLY by reading the LIVE consumer
(`_apply_walk` seeds from first-class ψ½ STATE post-route-(a), not
φ=0). Recipe: grep the retired NAMES across the WHOLE tree; refs tagged
"retired" are fine tombstones (L7); untagged live-framed refs — even in
`#` comments or non-automodule'd docstrings that produce NO build
warning — are Cardinal-Rule-1 staleness. Fix in-scope, FLAG out-of-scope
(this run left dangling live-framed `CarlsonInwardSweep`/
`AngularEdgeExtrapolation`/`CarlsonSweepContext` refs in
`augmented_mesh.py`, `solver.py`, `reduced_operator.py` for the main
agent — outside the doc-rider's file scope).

Complements [[feedback-double-category-architecture-insight]] (both
document an already-shipped type/algebra insight as a
success-resolution page) and the AGENT.md close-out arc (this is the
SUCCESS analogue — the deferral ledger is the success-story sibling of
the close-out's "infrastructure retained").
