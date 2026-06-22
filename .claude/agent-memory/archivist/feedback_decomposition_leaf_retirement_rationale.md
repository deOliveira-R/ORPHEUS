---
name: decomposition-leaf-retirement design-rationale expansion
description: Expanding a CONDENSED deep-dive (a `.. todo:: Archivist` marker) when #238 retired a DECOMPOSITION-into-separately-applicable-typed-leaves (StreamingOperator M_spatial / M_angular_redist over _SpatialSweepDirection / _MSpatialOperatorSum / AngularRedistributionOperator) into ONE fused in-sweep (L+C)ψ matvec — the survivor is the PHYSICS (WDD spatial + M-M angular recurrences, both sequentially coupled = MA-Q1 forbids TP) PLUS an honest "what was tried and retired" rationale; the orphan-smell-one-level-down argument; the five-anticipated-consumers list-table
type: feedback
---

A second-generation sibling of [[feedback_type_retirement_concept_survives]]
(S6.4(f)) and [[feedback_sweep_walk_collapse_relayering]] (S6.4(e)). The
distinction: those retired a TYPED FIELD or a method-product. #238
(`refactor/sn-foundation-cleanup`) retired a **DECOMPOSITION** — the
streaming matvec's two separately-applicable typed leaves (`M_spatial =
_SpatialSweepDirection(+1) + _SpatialSweepDirection(-1)` wrapped in
`_MSpatialOperatorSum`, plus the curvilinear `M_angular_redist` /
`AngularRedistributionOperator`) — folded into ONE fused in-sweep `(L+C)ψ`
matvec (`_OneDimScanWalk._apply_walk` in `loss_representation.py`). The
method-implementer had ALREADY condensed the three `wave-t-*-deep-dive`
sections and left a `.. todo:: Archivist` at `wave-t-orchestrated-apply`.
My job: EXPAND to full-narrative standard, NOT rebuild the architecture.

**The task shape was a STUB-EXPANSION (the todo marker), not a fresh
retirement sweep** — the de-role / past-tense / succession-note work was
mostly DONE by the method-implementer. Read the condensed version FIRST;
my value-add was the missing DERIVATION DEPTH and the HONEST RATIONALE,
not re-litigating the retirement.

**The two load-bearing expansions (Cardinal Rule 3):**

1. **The surviving PHYSICS, derived not asserted.** The condensed version
   stated the WDD recurrence and the M-M recurrence as bare equations. I
   added: (a) WHERE the WDD recurrence comes from (the 2-unknown cell
   balance + the DD mean-of-faces closure → solve for cell average); (b)
   WHY it forbids a TP (unroll to the lower-triangular `∏_{j<i}(...)`
   closed form — the off-diagonal products carry info BETWEEN cells AND
   depend on the ordinate μ through the denominators → spatial factor and
   angular index ENTANGLED → disjoint-axes contract fails); (c) the M-M
   angular term's geometric ORIGIN (the curved-trajectory frame-rotation
   derivative `(1-μ²)/r ∂ψ/∂μ`, cross-ref `:doc:discrete_ordinates` Eq.
   `balance-general` — DON'T duplicate the canonical M-M discretization,
   cross-ref it); (d) the ANGULAR-analogue framing: M-M half-grid
   recurrence is the angular twin of the spatial WDD obstruction, and the
   two run NESTED (spatial outer, angular inner) inside ONE walk — that
   nesting is precisely what SOTP cannot express (a flat sum of per-axis
   factors has no notion of one axis's recurrence running INSIDE
   another's). This nesting framing is the sharpest "why not a TP" content.

2. **The honest "design rationale / what was tried and retired"** — a
   FIVE-ROW `.. list-table::` "anticipated consumer | why a separate leaf
   was expected to help | what actually shipped (why it didn't
   materialise)". The five = Wave-O adjoint / #2 DSA / #200
   block-inverse-precond / per-direction debugging / slow per-direction
   test fallback. The KEY per-row arguments (these are the durable
   intellectual content, reuse verbatim-ish):
   - **Adjoint**: #240 G-adjoint rides the FUSED `loss_action_transpose`;
     transpose of a lower-triangular sweep is an upper-triangular sweep
     over the SAME coupling — splitting per-direction first gains nothing.
   - **DSA (#2)**: consumes the fused residual + a SEPARATE in-algebra
     diffusion op `A_diff = L+C−S`; never needs streaming-split-from-
     collision, let alone forward-split-from-backward.
   - **#200 block-inverse**: the full M-M sweep is ALREADY the natural
     O(N) EXACT inverse `(L+C)⁻¹` (one pass) → a block split would be
     WEAKER-NOT-CHEAPER (a block precond APPROXIMATES an inverse the sweep
     computes EXACTLY). This is the single most persuasive row.
   - **debugging / test-fallback**: a standalone direction summand that
     re-runs the whole bidirectional sweep + masks the other direction is
     a WORSE surface; the cross-check that mattered is the
     structurally-independent MMS, NOT a self-referential bit-identity
     between the fused pass and the sum of its own summands.

**The "orphan-smell one level down" argument is the close-out's spine.**
Frame it exactly: after #206 moved the walk into the loss-rep, the leaves
had ZERO production consumers and existed ONLY to be applied by their own
structural tests (the algebra-decomposition invariant `(L+C)ψ ≡
M_spatial·ψ + M_angular_redist·ψ`). A self-referential test loop —
machinery kept alive solely to feed the tests that exist solely to verify
that machinery — is a Cardinal-Rule-2 smell. Then the DECISIVE structural
observation (NOT incidental): the split could never be cheaper (the
forward-backward coupling forces a standalone per-direction apply to run
the whole bidirectional sweep anyway; the orchestrator existed only to
claw the 1.5× naive-sum penalty back to the fused-pass floor) NOR more
modular (exposing the directions as leaves forced them to SHARE NAMED
STATE — the `_compute_LpC` lifted outflow — i.e. a "modular" surface that
cannot be evaluated independently = the Pattern-4 illegal state).

**The "if a future consumer genuinely needs a leaf" `.. note::`** — close
the door against re-introduction: do NOT resurrect `_MSpatialOperatorSum`;
the structural obstruction (forward-backward coupling makes the directions
inseparable) must be addressed FIRST — the consumer needs a formulation
where the coupling ITSELF is the preconditioned object (in-algebra
diffusion op for DSA), not a re-split into directions that secretly share
state. Build against the fused `loss_action` / `loss_action_transpose`.

**Verification ground = the surviving MMS.** The retired
algebra-decomposition invariant only verified the split RECONSTRUCTED
ITSELF (a tautology, NOT a correctness claim). The surviving guarantee
that the curvilinear angular-redistribution term is correct is the
anisotropic curvilinear MMS (`:ref:sn-mms-curvilinear-aniso-verification`,
`catches("ERR-026")`, `tests/sn/verification/mms/test_curvilinear_aniso_convergence.py`)
— a structurally-independent L1 ground. SAY THIS explicitly in the
test-fallback row AND the subtraction-smell subsection.

**Preserved-anchor / preserved-label discipline (same as the family):**
- The 4 eq-labels (`wdd-forward-recurrence`, `mm-half-grid-recurrence`,
  `wave-t-cell-balance-three-terms`, `wave-t-mspat-curvilinear-subtraction`)
  carry `.. vv-status: <label> documented` and feed the auto-gen
  `docs/verification/matrix.rst` — KEEP verbatim, edit only surrounding
  prose. They land in matrix's Documented-only section (correct — they are
  structural/representational identities, not solver claims).
- The 3 section anchors — `wave-t-orchestrated-apply` is referenced by
  `docs/api/numerics.rst:182` (`:ref:`), so KEEP the anchor even though I
  RETITLED the section ("One bidirectional pass — design rationale and the
  retired per-direction split"). Retitle freely, keep anchor → cross-doc
  `:ref:` auto-picks the new title.

**Build / matrix delta:** MAIN-checkout `-E` baseline = 1 warning (mesh.py
`:paramref:` ERROR, out of scope; matches AGENT.md, NOT the worktree's 11).
Verification `-E` build EXIT=0, severity SET byte-identical, ZERO targets.
matrix.rst AUTO-REGENERATED by the build (do NOT hand-edit, do NOT stage
_build): `test_streaming_operator` foundation 67→54 (−13, the deleted
decomposition-class tests), Total 5358→5345. Brief estimated "5 fewer";
actual was 13 — report the REAL number, it's still "the deleted
decomposition classes, expected". The 4 preserved labels STILL indexed
in the regenerated matrix.

**Em-dash title-underline pitfall (recurs every time):** new H3 subsection
titles with em-dashes (`—`, 1 code point, 3 bytes) — size underline with
`len(title)` in python NOT `wc -c`. I ran a python pass to fix; it ALSO
"fixed" 10 PRE-EXISTING tolerated over-run underlines elsewhere in the
file (RST allows underline LONGER than title with no warning) → those were
OUTSIDE my scope, so I RESTORED them (kept the diff to only my 3 new
underlines). Lesson: a global underline-normalize pass touches unrelated
lines; scope it to the lines you added, or restore the rest.

**File marker ladder:** `operator_algebra.rst` uses `=`/`-`/`~`/`^`; `~`
IS an established H3 under `----` H2 (line 2620 etc.) → my new `~~~`
subsections are the correct level. Scan first-appearance markers before
picking.
