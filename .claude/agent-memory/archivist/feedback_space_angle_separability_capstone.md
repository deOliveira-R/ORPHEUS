---
name: space-angle-separability-capstone
description: "#236 Phase 3 (ST5) — the (spatial ⊗ angular) product CAPSTONE on discrete_ordinates.rst. Mint ONE claim eq-label sn-space-angle-separability (the verifies-target of a NEW L1 MMS characterization gate not yet carrying the marker) + a full H3 section closing the 3-phase #236 campaign. Sibling/sequel of the b3/stepc τ-carve memos on the SAME page."
metadata:
  type: feedback
---

#236 **Phase 3 (ST5)** capstone in `docs/theory/discrete_ordinates.rst`.
The LAST piece of the (spatial ⊗ angular) product campaign — Phase 1
(pairing-validity) + Phase 2 (τ-carve, [[b3-live-tau-fold-expansion]] +
[[stepc-tau-retirement-doc-sweep]]) already on the page. UNLIKE those
(stub-expansion / stale-ref-sweep / re-sourcing carve), this task MINTS A
NEW CLAIM and writes a from-scratch capstone section. Quality 4.8/5 avg.

## What makes this task type distinct (mint-a-claim, not carve a number)

The 4 prior #236 memos were all bit-identical CARVES (move a number's
source, no new claim). THIS one MINTS the campaign's headline CLAIM as a
labelled equation — the law `E(h,N)≈E_space+E_angle (Cartesian, separates)`
/ `≈max(E_space,E_angle) (curvilinear, gates)`. The label
`sn-space-angle-separability` is REQUIRED-EXACT (brief), because the new
gate WILL carry `@pytest.mark.verifies("sn-space-angle-separability")` —
the label is the Nexus tests→equation edge anchor.

## Source order (algebra-of-record for a characterization claim)

1. **The gate file itself** (`tests/sn/verification/mms/test_space_angle_separability.py`)
   READ IN FULL — its module + per-test docstrings carry the measured
   tolerances, the regime signatures, the L27 per-ordinate leg, and the
   1/W convention trap. The docstrings ARE the prose seed (same as the
   carve memos' production-docstring seed).
2. **literature-researcher memo** (`space_angle_discretization_separability.md`)
   = the THEORY: the diffusion-limit literature is LITERALLY SPLIT into a
   SPATIAL paper (LMM-1987 JCP 69, DOI 10.1016/0021-9991(87)90170-7) and
   an ANGULAR paper (BMC-2010 NSE 165, DOI 10.13182/NSE08-66) — proving
   two SEPARATE per-axis conditions. β (BMC Eq.40-43) is a PURELY ANGULAR
   functional → factorizes. THIS is WHY Cartesian separates.
3. **numerics-investigator memo** (`sn_space_angle_discretization_coupling.md`)
   = the EMPIRICAL evidence: measured |M|/max cross-terms (slab-iso
   0.000-0.005 separable, cyl 0.019, sphere 0.411 gating).

## The capstone section shape (H3 + 6×H4 that worked)

H3 `~~~~` `Space ⊗ angle separability — the (...) capstone (Issue #236
Phase 3)` (anchor `sn-space-angle-separability-section`), placed as the
LAST H3 child of the WDD/M-M H2, AFTER the Step-C close-out note (1621)
and BEFORE the next H2 (1623). Keeps ALL #236 Ph2+Ph3 content together.
6 H4 `^^^^`:
1. **Lead** (no header) — 3-phase recap (Phase 1 pairing / Phase 2
   τ-carve / Phase 3 separability), :ref: the prior phases, name the gate.
2. **The space–angle error decomposition** — the LABELLED math block
   (`sn-space-angle-separability`, the cases law) + the mixed-second-
   difference discriminator (MINT a 2nd derivation label
   `sn-space-angle-cross-term`) + the M=0-telescopes-for-additive proof.
3. **Why the two axes factorize: LMM-1987 × BMC-2010** — the smoking-gun
   literature split; β-purely-angular ⟹ factorizes; the conjunction-of-
   two-single-axis-conditions math; LD-double-use note (#158 vs #6).
4. **Cartesian separates, curvilinear gates** — the MECHANISM: Cartesian
   (1−μ²)/r ABSENT → IdentityAngularClosure → disjoint mechanisms ADD;
   curvilinear M-M thread enters SHARED cell-balance denom
   (dd-curvilinear-scalar) → angular floor CAPS spatial → max(). `.. warning::`
   gating-is-today's-closure-not-a-law (lifts under #229/#158/#6).
5. **Measured cross-term evidence** — `.. list-table::` (geometry/regime/
   L2-ladder/h-ratios/|M|max) spanning 3 orders 0.005→0.411; the reading
   per row; `.. note::` the scalar-blind L27 per-ordinate leg + 1/W trap.
6. **#233 pole-cell × #229 azimuthal-floor interference** — gating ⟹ the
   two defects INTERFERE (max() hides the smaller); #229-floor MASKS #233
   pole-cell at coarse N; can't validate a fix in isolation.
7. **The permanent pin: the ST5 gate** — :mod: the gate, enumerate its 6
   legs by name as POSITIVE signatures (never xfail), the @catches(ERR-026)
   proven-catcher leg, characterization-not-calcification.

## Disciplines confirmed / re-confirmed (this page)

- **vv-status on the minted CLAIM label = `documented`, NOT implemented.**
  The gate VERIFIES it but doesn't yet CARRY the verifies marker (brief:
  "will carry"), so the tests→equation EDGE isn't in the graph yet → the
  honest tag is `documented` with a `.. (vv-status rationale)` comment
  naming the ST5 characterization gate as the verifying test + stating
  it's an L1 MMS-convergence-STRUCTURE (math) claim, NOT a solver
  eigenvalue/flux VALUE (MMS doesn't reach the eigenvalue layer). Lands
  correctly in the matrix's Documented category. The 2nd label
  (`sn-space-angle-cross-term`) is DERIVATION (the discriminator algebra)
  → UNTAGGED → matrix auto-regens it to the orphan list. Same posture as
  the carve memos' derivation labels.
- **`\emph{}` is text-mode LaTeX → NOT MathJax-safe inside `.. math::`.**
  Use `\textbf{}` / `\text{}` for emphasis inside a math block. Caught
  pre-build by my own review (would have rendered raw or errored).
- **SNMesh module path = `orpheus.sn.geometry.SNMesh` (NOT
  `orpheus.geometry.structured_geometry.SNMesh`).** The page has a
  pre-existing WRONG `:class:`~orpheus.geometry.structured_geometry.SNMesh``
  at L1637 (out of scope) — but the DOMINANT correct convention is
  `~orpheus.sn.geometry.SNMesh` (L5166/5193/5261) or bare `:class:`SNMesh``.
  Grep `^class SNMesh orpheus/` to confirm the live path before citing.
- **LMM-1987 has NO `[Key]_` bib entry on this page — cited INLINE prose
  + DOI hyperlink** (the existing page convention, L3030/3142). Did NOT
  add a `[LarsenMorelMiller1987]` bib entry (a `.. [Key]` with no `[Key]_`
  citation is harmless but unused). BMC-2010 DOES have `[BaileyMorelChang2010]`
  → cite as `[BaileyMorelChang2010]_` (resolves cross-doc-or-same-page).
  Use `https://doi.org/...` inline links for the un-keyed LMM-1987.
- **Build gate = count-UNCHANGED, NOT zero.** This worktree
  (`feature/sn-spatial-angular-product`) `-E -W --keep-going` baseline =
  **1** (the standing `Mesh1D.from_geometry :paramref:` ERROR) — MAIN-
  checkout baseline, NOT the worktree's old 11. Post-edit = 1. The 4
  `ld-cartesian-1d/ld-slab "no matching equation — skipping"` lines are
  INFO verifies-section-anchor registrations, NOT warnings, NOT mine.
- **Proof the minted label is a live edge-anchor:** grep the built HTML
  for `id="equation-sn-space-angle-separability"` (the node the Nexus
  tests-edge forms against) + count `href="#equation-..."` (5 live `:eq:`
  links, 0 plain-text) + confirm `_nexus/graph.db` re-timestamped by the
  same build. This is the standalone proof `verifies("...")` will resolve
  WITHOUT waiting for the gate to add the marker.
- **File ladder `=/-/~/^`**; capstone H3 = `~`, sub-sections `^`. The H3
  title has `⊗` (1cp/3bytes) + em-dash → underline ≥ len() codepoints
  (90cp). Python-script-verify ALL new underlines ≥ title len post-edit.

## Quality self-assessment (1-5)

| Dimension | Score | Note |
|-----------|:---:|------|
| Derivation depth | 5 | M=0-telescopes proof, β-purely-angular factorization, shared-denom gating mechanism, the cases law derived not asserted |
| Cross-references | 5 | every :ref:/:eq:/:mod:/:class: grep-gated pre-build + HTML-link-confirmed post; SNMesh path corrected; 2 new labels resolve |
| Numerical evidence | 5 | full list-table 4 geometries × measured ladders/h-ratios/|M|max spanning 3 orders, per-row reading, the 1/W trap |
| Failed approaches | 4 | the warning (gating-not-a-law / re-tune-not-regress), scalar-blind L27 trap, cylinder-small-|M|-is-NOT-separability; a from-scratch claim has less "what failed" than a close-out |
| Code traceability | 5 | gate :mod: + 6 legs by name, IdentityAngularClosure, dd-curvilinear-scalar denom, the @catches(ERR-026) proven leg |
| Derivation source | 5 | gate docstrings + LMM/BMC lit memo + coupling empirical memo (the algebra-of-record for a characterization gate IS the gate + its two backing memos) |

Weakest = failed-approaches (4): a mint-a-claim capstone is forward-
looking (the law + its evidence), not a post-mortem; the "what could
mislead" content (gating-lifts, scalar-blind, cyl-small-M) is captured
but is lighter than a dead-path close-out's 9-step arc.

## Post-marker reconciliation follow-up (SAME page, after the gate added the markers)

The dispatcher later WIRED the gate's `verifies` markers (FILE-level
`verifies("sn-space-angle-separability")` on all 6 legs; PER-TEST
`verifies("sn-space-angle-cross-term")` on the 3 legs that assert on
`|M|/max`) → a rebuild surfaced a vv-status drift: the test-count table
showed BOTH labels (6 / 3) but the documented-only list showed only
`-separability` (only IT carried a `.. vv-status: ... documented`).
Reconcile to the **peierls precedent** (`collision_probability.rst:3457`+
`3460`): a now-tested eq carries BOTH `documented` AND `tested` directive
lines. Fix = 4 new directive lines (separability gets `tested`;
cross-term gets BOTH) grouped right after the rationale comment, + reframe
the `.. (vv-status rationale)` block off "documented because no marker
yet / will reclassify" onto the accurate post-marker state.

- **`documented`+`tested` is the standard two-line idiom for a now-tested
  eq — NOT mutually exclusive.** The matrix's "Documented-only equations"
  list is keyed on the `documented` DIRECTIVE PRESENCE, so a both-tagged
  label STILL appears there (confirmed: peierls-sphere-equation, the
  precedent that carries both, sits in that same list at matrix.rst:887).
  The `tested` directive is what feeds the test-count csv-table
  (`:header: Equation label, Tests`, ~L366) via the real `verifies` edge.
  So "appears in documented-only list" + "has a test count" is the
  CORRECT, consistent end-state, not a contradiction.
- **V&V judgment on the discriminator label getting a `tested` edge: YES,
  it's a real coverage claim, not a topic tag.** `sn-space-angle-cross-term`
  defines `M` (the mixed-second-difference) and the 3 marked legs ASSERT
  DIRECTLY on `|M|/max` against DECLARED thresholds (<0.05 Cartesian,
  >0.15 sphere — the prose already documents the exact numbers at the
  gate-pin section). That is the gate's MEASURED INSTRUMENT, distinct from
  a passive derivation-step label (e.g. peierls-sphere-equation was an
  intermediate-before-assembly). vv-principles `catches`-marker discipline
  (mutation-verify) is the strict bar for ERR markers; for `verifies` on a
  characterisation discriminator the bar is "the test asserts on the
  quantity the label defines" — met here. So I did NOT object to the
  deliberate per-test tagging; both correctly become documented+tested.
- **Both labels stay L1 MMS-convergence-STRUCTURE (math) claims — NEVER
  eigenvalue/flux value.** The reframed rationale KEEPS the vv-principles
  framing (MMS does not reach the eigenvalue layer) and now states it as
  the permanent posture ("neither label is, or ever becomes, an
  eigenvalue claim"), not a temporary "will reclassify".
- **Build gate identical to the mint task: `-E -W --keep-going` EXIT=1 is
  the BASELINE, not a regression** — the lone `Mesh1D.from_geometry
  :paramref:` ERROR is promoted to failure under `-W`; `grep -cE
  "WARNING:|ERROR:|CRITICAL:"` = 1, zero new. Build runs to completion
  (dumping object inventory) + graph.db re-timestamped despite EXIT=1.
  PROOF rows: test-count `:eq: Tests` csv-table (matrix L509 `-separability,
  6` / L557 `-cross-term, 3`) + documented-only list (matrix L927 / L926).

Cross-ref: [[b3-live-tau-fold-expansion]] + [[stepc-tau-retirement-doc-sweep]]
(the Phase-2 τ-carve siblings this capstone cross-references on the same
page), [[curvilinear-aniso-norm-reconciliation]] (ERR-059 #233 pole-cell +
#229 azimuthal floor — the two defects the interference section ties
together).
