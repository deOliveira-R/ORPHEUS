---
name: issue-138-step5b-first-consumer-closeloop
description: Doc-pass recipe for a "latent-consumer-now-realized" close-the-loop when a landed follow-on wires the FIRST production consumer of a verified-but-unwired operator type (+ the spectral-invisibility object-gate V&V section)
metadata:
  type: feedback
---

The #226 taxonomy-step doc-pattern family, step 5b (task #138): when a
landed follow-on wires the **first production consumer** of an operator
type that shipped "verified but not yet wired" (here `MatrixInverseOperator`,
first consumed by `solve_homogeneous_infinite`), the doc pass is a
**close-the-loop** across THREE surfaces. Instance: `refactor/inverse-as-operator`,
homogeneous.rst + operator_algebra.rst + api/numerics.rst.

**Why:** the follow-on's motivating doc (the operator's "consumer ruling"
section) literally said "the full spelling `K = MatrixInverseOperator(loss) @ F`
is the follow-on (task #138), **not this step**". Closing that loop is the
L-007 deferred-seam-LANDED tense-flip, applied across the whole blast radius.

**How to apply:**

- **Three-page tense-flip, each with a different job.** (1) The THEORY page
  of the consumer (homogeneous.rst) = the main rewrite: re-spell the eigensolve
  in the operator algebra (`K = MatrixInverseOperator(loss) @ production`;
  `dominant_eigenpair(K.as_matrix())`). (2) The OPERATOR's own consumer-ruling
  section (operator_algebra.rst) = flip "the follow-on (task #138), not this
  step" → "landed, first production consumer". (3) The API narrative
  (api/numerics.rst) = fix the stale "the solver uses `direct_eigenvalue`"
  claim. The theory-page grep is the FLOOR; the parallel API + operator pages
  go stale identically (L-001 whole-`docs/`-tree blast radius).

- **Don't over-flip "not yet wired" → "wired". Distinguish the AXES.** The type
  is now wired as a *production spelling* (explicit construction) — but the
  *factory routing* (auto `.inverse()` returning it) STILL waits, and the
  consumer DELIBERATELY bypasses it. Keep that distinction: "wired as a
  spelling: yes / wired via factory: no, and permanently bypassed here".

- **Strategy-choice-as-type is the load-bearing design point.** Explicit
  `MatrixInverseOperator(loss)` vs structure-keyed `loss.inverse()` (which,
  reading the sum `C − K_iso` with invertible leading term, returns the
  ITERATIVE `GreenOperator` splitting). The direct-realization decision encoded
  as a TYPE, not a `strategy=` flag — the §3 strategy-override seam realized
  honestly. Document WHY the iterative splitting is wrong for a 0-D dense block.

- **Shared-extraction-home reframe: sibling-spelling, not retirement.** When
  the follow-on extracts a shared primitive (`dominant_eigenpair` = the PF
  eig-extraction: argmax-real + complex-reject + sign-normalize), the old
  convenience engine (`direct_eigenvalue`) becomes the `(A,F)`-posed SIBLING
  that ALSO delegates to the shared home. Reframe it as sibling-spelling; it
  keeps **zero production consumers** as a non-superseded sibling engine + RQI
  test oracle (explicitly NOT the fuller-view-oracle exception — a peer engine,
  not a relinquished fuller view). VERIFY zero-consumers against live code
  (`grep 'direct_eigenvalue('` in orpheus/ minus the def).

- **The spectral-invisibility → object-level gate V&V section (a reusable
  teaching unit for ANY resolvent/eigenvalue doc).** Two mutations are
  *spectrally invisible* (|Δk| = 0 exactly, every k-gate blind): factor-order
  swap `F·A⁻¹` (SIMILAR to `A⁻¹F` via `A·M·A⁻¹ = F·A⁻¹`, similar ⇒ same
  spectrum; the eigenvector remaps `φ ↦ Aφ`) and resolvent TRANSPOSE
  (`eig(Mᵀ)=eig(M)` via char-poly transpose invariance). The committed catcher
  is the OBJECT-level matrix-identity gate (`K.as_matrix ≡ np.linalg.solve(A,F)`,
  rtol=1e-12), where both move `[K]` by O(1). Lesson: **pin the object, not
  just its spectrum** — a value gate can be blind to a whole mutation class for
  structural (spectral-similarity) reasons. Give the similarity identity an
  eq-label + `.. vv-status: <label> documented` (structural, not a solver
  claim → lands in matrix Documented-only). Candidate vv-principles uplift.

- **Verification mechanics (all clean here):** matrix auto-regens (L-008) —
  `numerics/test_eigenvalue` 30→39, `homogeneous/test_homogeneous` 12→14,
  `removal-matrix` verifiers 61→63 — NOT warnings. `dominant_eigenpair` `:func:`
  refs render as LINKS (eigenvalue.py IS `automodule`'d at api/numerics.rst:480,
  0 `:label:` docstrings → no collision). Preserve every verifies-target label
  (removal-matrix / fission-matrix / fixed-source-solve / keff-update /
  matrix-eigenvalue): the equation BODIES stay (`M=A⁻¹F` unchanged), only the
  surrounding prose flips (L-003). Build `-E -W` EXIT 0, WARNING/ERROR/CRITICAL
  set unchanged (0).
