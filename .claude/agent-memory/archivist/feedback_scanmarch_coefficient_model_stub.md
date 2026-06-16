---
name: scanmarch-coefficient-model-stub-expansion
description: #240 Phase 2 Step D5a — expand loss-rep-scanmarch-coefficient-model stub + REWRITE two STALE verifies-target equations (scanmarch-solve/-apply) into coefficient-model form. SIBLING of affine-in-sigma (same loss_representations.rst page, same #240 Ph2). The stale eqs are TEST-PINNED labels → keep label, rewrite body only.
metadata:
  type: feedback
---

#240 Phase 2 Step D5a (commit `66dbd9a`): the 2-D ScanMarch row-march's inline
DiamondDifference folded onto the scheme coefficient model. Two-part docs task
on `docs/theory/loss_representations.rst`: (1) expand the
`loss-rep-scanmarch-coefficient-model` `.. todo::` stub; (2) rewrite two
now-STALE equations `loss-rep-scanmarch-solve` / `loss-rep-scanmarch-apply`.
SIBLING of [[affine-in-sigma-stub-expansion]] (SAME page, SAME #240 Phase 2).

**Why this one is distinct from the affine sibling — the STALE-equation rewrite
is the load-bearing half, and the two stale labels are TEST-PINNED.**

- **The closeout memo IS the algebra-of-record** here (not a test docstring like
  the affine sibling). `method-implementer/issue_240_phase2_step_d5a_closeout.md`
  carried the full fold mechanism + the LOAD-BEARING FINDING (the 2-D `_sweep_interior`
  was the ONE remaining `÷D` division path → re-baseline ~1 ULP; 1-D scan
  byte-identical = negative control). But STILL READ THE CODE @ `66dbd9a` to get
  exact coefficient algebra: `diamond.py:_cartesian_streaming_diagonal`
  (`couplings[a]=2g_a`, `denom=Σ_t+Σ2g_a`, explicit LEFT FOLD), `cartesian_scan_coefficients`
  (`a=2·scan_diag·inverse_denom−1`, `scan_diag=2g_x`), `reflect_scan_coefficients`
  (`α=−1,β=2ψ̄` via `_reflection_coeffs(ψ̄,_DD_W)`), `scheme.py:source_emission`
  (`β=QV·inverse_denom/w`), `cell_average`/`outgoing_face_from_average`;
  `loss_representation.py` `_sweep_interior`/`_loss_action_interior`;
  `scan.py:_scanmarch_row` (gained `w` param).

**The STALE-equation-rewrite discipline (the CRITICAL half):**

- **GREP THE MATRIX FIRST.** `docs/verification/matrix.rst` listed BOTH stale
  labels (`loss-rep-scanmarch-solve`, `-apply`) at "count 11" = they are
  `@pytest.mark.verifies(...)` TARGETS (one test:
  `test_scan_march_equivalence.py` verifies `loss-rep-scanmarch` + `-solve` +
  `-apply` = the ScanMarch≡FullFieldWavefront oracle). ⟹ **KEEP both label
  NAMES, rewrite only the equation BODY.** The claim (scan-march coefficients +
  their solve/apply realizations) is unchanged; only the inline-DD→coefficient-
  model presentation changed; the oracle still pins it. Renaming a verifies-
  target label would orphan the test's `verifies(...)`.
- **A stale-equation rewrite ripples to every `:eq:` CITING it.** The old eqs
  defined `s_a = 2|μ_a|/Δa` (DD-prescaled). The coefficient-model rewrite makes
  the RAW `g_a = |μ_a|/Δa` primary + the scheme-scaled coupling `c_a = 2g_a`.
  THREE other in-page sites cited the old `s_a`/`β=2(...)/D`: the affine-carve
  `loss-rep-affine-cell` block (`denom=Σ_a s_a + σ`), its sweep-contrast
  paragraph (`denom=Σ_a s_a + σ`), and the facewise-section prose (`s_y ψ_y`,
  `β=2(Q+s_y ψ_y)/D`, the DD-facewise `single 0th-order face value`). RENAMED
  all `s_a → c_a` for one-page consistency + retargeted the facewise prose's
  `:eq:` to the new `-solve-affine` sub-label. **A rename in a definitional eq
  is a whole-page sweep, not a local edit.**
- **SPLIT a busy stale eq into a label-preserving primary + a NEW sub-label.**
  Old `loss-rep-scanmarch-solve` crammed `α`, `β`, `D` into one block. Rewrote
  as: `-solve` (KEPT label) = the coupling/diagonal `c_a=2g_a, S=Σ_t+c_x+c_y`;
  NEW `-solve-affine` = the `(α,β)` scan coefficients. Same for apply: `-apply`
  (KEPT) = the reflection `(α=−1,β=2ψ̄)`; NEW `-apply-residual` = the residual.
  The KEPT labels stay the verifies-targets; the NEW labels are derivation
  decompositions (orphan-list, untagged — see below).

**New derivation eq-labels → orphan-list, UNTAGGED (sibling precedent).** The
two new labels (`-solve-affine`, `-apply-residual`) auto-appeared in
matrix.rst's "Orphan equations" list (65→67) right beside the affine sibling's
own untagged derivation labels (`loss-rep-affine`/`-affine-cell`/`-leaf-sum`/
`-removal-sigma`). This is CORRECT: derivation/representational decompositions
of an ALREADY-VERIFIED claim (the oracle pins `-solve`/`-apply`) are NOT
independent solver claims → leave UNTAGGED (do NOT `:vv-status: documented`-tag),
matching the four established siblings. matrix.rst is AUTO-GENERATED → the
65→67 bump is expected auto-regen, NEVER hand-edit. The orphan list is NOT a
build warning (build stayed clean).

**Cross-ref path accuracy (Cardinal Rule 1, even though all render plain-text):**
the whole page's `:meth:`/`:class:` render plain-text (page convention; not
member-autodoc'd), so path is invisible in OUTPUT but the RST SOURCE must name
the OWNING class. Runtime-checked: `cartesian_scan_coefficients`,
`reflect_scan_coefficients`, `residual_kernel_batch`, `cell_kernel_batch`,
`affine_scan_coefficients`, `source_emission`, `cell_average`,
`outgoing_face_from_average` ALL live on `DiscretizationSchemeBase` ONLY (the
Protocol `DiscretizationScheme` declares ONLY `update`/`residual`; the stub's
`DiscretizationScheme.cartesian_scan_coefficients` was WRONG). Convention used
(matches the committed page's `cell_kernel_batch` refs): DD-overridden
producers/kernels → `diamond.DiamondDifference.<m>`; generic reconstruction ops
(`source_emission`/`cell_average`/`outgoing_face_from_average`) →
`scheme.DiscretizationSchemeBase.<m>`. `_reflection_coeffs` is a DD staticmethod
→ `diamond.DiamondDifference._reflection_coeffs` (the stub had it on `scheme.`).

**Section-marker ceiling.** The stub's parent
(`loss-rep-scanmarch-coefficient-model`) is already H3 `~~~~`; `^^^^` (H4) is
NOT used ANYWHERE in this file (ladder = `=`/`-`/`~`). So the rich narrative
stayed SINGLE-LEVEL (one H3) with bold-lead paragraphs + admonitions + 3
`.. list-table::` (two-doors / apply-vs-solve-asymmetry / verification-gates) —
NO deeper headings (would trip "Inconsistent title style: skip level N→N+2").
This differs from the affine sibling which had `~~~~` sub-subsections under an
H2 `----` parent.

**Production docstring staleness flagged, OUT of scope.** `_sweep_interior`'s
docstring (`loss_representation.py:1302-1304`) STILL carries the pre-D5a inline
form `α = 2s_x/D − 1, β = 2(Q + s_y ψ_y)/D` in PROSE (the body is coefficient-
model-correct). Task said do NOT edit production code → flagged for main agent,
not touched.

**Build gate:** `-E -W --keep-going` baseline = EXIT 1 / 1 ERROR (the
pre-existing `orpheus/geometry/mesh.py` `Mesh1D.from_geometry` `paramref` role —
`sphinx-paramlinks` NOT installed/configured → adding the dep is an env decision,
OUT of scope for a docs task). Post-edit IDENTICAL (count-unchanged = pass). My
edits added ZERO new WARNING/ERROR/CRITICAL. Verified rendered HTML: all 4
eq-nodes present, section heading rendered, `.. todo::` gone.

**Quality self-assessment (1-5):**

| Dimension | Score | Note |
|-----------|:---:|------|
| Derivation depth | 5 | full coefficient-model fold: raw-g-in/scaled-c-out, scheme-door table, solve closes-generically-in-w, apply reflection-scan+residual-kernel, the division-to-reciprocal re-baseline mechanism |
| Cross-references | 4 | every :eq:/:ref: resolves + path-corrected to owning class; :meth:/:class: render plain-text by page convention (not a regression) |
| Numerical evidence | 5 | re-baseline admonition: oracle k_inf=1.875/φ=Q/Σ_t + test_keff_2d closed-form k_inf=νΣ_f/Σ_a + SI≡Krylov 2.8e-12 + 1-D byte-identical negative control; verification-gate table |
| Failed approaches | 4 | the inline-DD that worked-for-DD-only, the ÷D-was-the-last-division finding, #242 affine-merge deferral, honest-interim D5a/D5b split — but this is a SUCCESS fold not a falsification, so less failed-history than a close-out |
| Code traceability | 5 | every method named to its owning class, _DD_W/_cartesian_streaming_diagonal/_reflection_coeffs single-sourcing, every gate cited |
| Derivation source | 5 | closeout memo (algebra-of-record) cross-checked against code @ 66dbd9a; equations match _cartesian_streaming_diagonal/cartesian_scan_coefficients/reflect_scan_coefficients/residual_kernel_batch exactly |

Weakest = cross-references (4) only because of the page-wide plain-text code-xref
convention (surfacing the whole scheme/operator package is its own arch-docs task).
